# coding=utf-8

import pytest

from datetime import datetime, timedelta

import numpy as np
import pandas as pd

from pandas import Series, NaT, date_range, Timestamp, Timedelta

from pandas.compat import range
import pandas.util.testing as tm


@pytest.fixture
def tdser():
    return Series(['59 Days', '59 Days', 'NaT'], dtype='timedelta64[ns]')


class TestTimedeltaSeriesArithmeticWithIntegers(object):
    # Tests for Series with dtype 'timedelta64[ns]' arithmetic operations
    # with integer and int-like others

    # ------------------------------------------------------------------
    # Addition and Subtraction

    def test_td64series_add_int_series_invalid(self, tdser):
        with pytest.raises(TypeError):
            tdser + Series([2, 3, 4])

    @pytest.mark.xfail(reason='GH#19123 integer interpreted as nanoseconds')
    def test_td64series_radd_int_series_invalid(self, tdser):
        with pytest.raises(TypeError):
            Series([2, 3, 4]) + tdser

    def test_td64series_sub_int_series_invalid(self, tdser):
        with pytest.raises(TypeError):
            tdser - Series([2, 3, 4])

    @pytest.mark.xfail(reason='GH#19123 integer interpreted as nanoseconds')
    def test_td64series_rsub_int_series_invalid(self, tdser):
        with pytest.raises(TypeError):
            Series([2, 3, 4]) - tdser

    @pytest.mark.parametrize('scalar', [1, 1.5, np.array(2)])
    def test_td64series_add_sub_numeric_scalar_invalid(self, scalar, tdser):
        with pytest.raises(TypeError):
            tdser + scalar
        with pytest.raises(TypeError):
            scalar + tdser
        with pytest.raises(TypeError):
            tdser - scalar
        with pytest.raises(TypeError):
            scalar - tdser

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [
        np.array([1, 2, 3]),
        pd.Index([1, 2, 3]),
        pytest.param(Series([1, 2, 3]),
                     marks=pytest.mark.xfail(reason='GH#19123 integer '
                                                    'interpreted as nanos'))
    ])
    def test_td64series_add_sub_numeric_array_invalid(self, vector,
                                                      dtype, tdser):
        vector = vector.astype(dtype)
        with pytest.raises(TypeError):
            tdser + vector
        with pytest.raises(TypeError):
            vector + tdser
        with pytest.raises(TypeError):
            tdser - vector
        with pytest.raises(TypeError):
            vector - tdser

    # ------------------------------------------------------------------
    # Multiplicaton and Division

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [np.array([20, 30, 40]),
                                        pd.Index([20, 30, 40]),
                                        Series([20, 30, 40])])
    def test_td64series_div_numeric_array(self, vector, dtype, tdser):
        # GH 4521
        # divide/multiply by integers
        vector = vector.astype(dtype)
        expected = Series(['2.95D', '1D 23H 12m', 'NaT'],
                          dtype='timedelta64[ns]')

        result = tdser / vector
        tm.assert_series_equal(result, expected)

        with pytest.raises(TypeError):
            vector / tdser

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [np.array([20, 30, 40]),
                                        pd.Index([20, 30, 40]),
                                        Series([20, 30, 40])])
    def test_td64series_mul_numeric_array(self, vector, dtype, tdser):
        # GH 4521
        # divide/multiply by integers
        vector = vector.astype(dtype)

        expected = Series(['1180 Days', '1770 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        result = tdser * vector
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [
        np.array([20, 30, 40]),
        pytest.param(pd.Index([20, 30, 40]),
                     marks=pytest.mark.xfail(reason='__mul__ raises '
                                                    'instead of returning '
                                                    'NotImplemented')),
        Series([20, 30, 40])
    ])
    def test_td64series_rmul_numeric_array(self, vector, dtype, tdser):
        # GH 4521
        # divide/multiply by integers
        vector = vector.astype(dtype)

        expected = Series(['1180 Days', '1770 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        result = vector * tdser
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('one', [1, np.array(1), 1.0, np.array(1.0)])
    def test_td64series_mul_numeric_scalar(self, one, tdser):
        # GH 4521
        # divide/multiply by integers
        expected = Series(['-59 Days', '-59 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        result = tdser * (-one)
        tm.assert_series_equal(result, expected)
        result = (-one) * tdser
        tm.assert_series_equal(result, expected)

        expected = Series(['118 Days', '118 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        result = tdser * (2 * one)
        tm.assert_series_equal(result, expected)
        result = (2 * one) * tdser
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('two', [
        2, 2.0,
        pytest.param(np.array(2),
                     marks=pytest.mark.xfail(reason='GH#19011 is_list_like '
                                                    'incorrectly True.')),
        pytest.param(np.array(2.0),
                     marks=pytest.mark.xfail(reason='GH#19011 is_list_like '
                                                    'incorrectly True.')),
    ])
    def test_td64series_div_numeric_scalar(self, two, tdser):
        # GH 4521
        # divide/multiply by integers
        expected = Series(['29.5D', '29.5D', 'NaT'], dtype='timedelta64[ns]')

        result = tdser / two
        tm.assert_series_equal(result, expected)


class TestTimedeltaSeriesArithmetic(object):
    def test_td64series_add_sub_timestamp(self):
        # GH11925
        tdser = Series(pd.timedelta_range('1 day', periods=3))
        ts = Timestamp('2012-01-01')
        expected = Series(date_range('2012-01-02', periods=3))
        tm.assert_series_equal(ts + tdser, expected)
        tm.assert_series_equal(tdser + ts, expected)

        expected2 = Series(date_range('2011-12-31', periods=3, freq='-1D'))
        tm.assert_series_equal(ts - tdser, expected2)
        tm.assert_series_equal(ts + (-tdser), expected2)

        with pytest.raises(TypeError):
            tdser - ts

    def test_timedelta64_operations_with_DateOffset(self):
        # GH 10699
        td = Series([timedelta(minutes=5, seconds=3)] * 3)
        result = td + pd.offsets.Minute(1)
        expected = Series([timedelta(minutes=6, seconds=3)] * 3)
        tm.assert_series_equal(result, expected)

        result = td - pd.offsets.Minute(1)
        expected = Series([timedelta(minutes=4, seconds=3)] * 3)
        tm.assert_series_equal(result, expected)

        result = td + Series([pd.offsets.Minute(1), pd.offsets.Second(3),
                              pd.offsets.Hour(2)])
        expected = Series([timedelta(minutes=6, seconds=3), timedelta(
            minutes=5, seconds=6), timedelta(hours=2, minutes=5, seconds=3)])
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

    def test_operators_timedelta64(self):
        # series ops
        v1 = date_range('2012-1-1', periods=3, freq='D')
        v2 = date_range('2012-1-2', periods=3, freq='D')
        rs = Series(v2) - Series(v1)
        xp = Series(1e9 * 3600 * 24,
                    rs.index).astype('int64').astype('timedelta64[ns]')
        tm.assert_series_equal(rs, xp)
        assert rs.dtype == 'timedelta64[ns]'

        df = pd.DataFrame(dict(A=v1))
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

    def test_td64_sub_NaT(self):
        # GH#18808
        ser = Series([NaT, Timedelta('1s')])
        res = ser - NaT
        expected = Series([NaT, NaT], dtype='timedelta64[ns]')
        tm.assert_series_equal(res, expected)

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

    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_operators_timedelta64_with_timedelta_invalid(self, scalar_td):
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # check that we are getting a TypeError
        # with 'operate' (from core/ops.py) for the ops that are not
        # defined
        pattern = 'operate|unsupported|cannot'
        with tm.assert_raises_regex(TypeError, pattern):
            td1 * scalar_td
        with tm.assert_raises_regex(TypeError, pattern):
            scalar_td * td1
        with tm.assert_raises_regex(TypeError, pattern):
            scalar_td ** td1
        with tm.assert_raises_regex(TypeError, pattern):
            td1 ** scalar_td

    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_timedelta_rfloordiv(self, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan
        result = scalar_td // td1
        expected = Series([1, 1, np.nan])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_timedelta_rfloordiv_explicit(self, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # We can test __rfloordiv__ using this syntax,
        # see `test_timedelta_rfloordiv`
        result = td1.__rfloordiv__(scalar_td)
        expected = Series([1, 1, np.nan])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_timedelta_floordiv(self, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        result = td1 // scalar_td
        expected = Series([0, 0, np.nan])
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_td64_series_with_tdi(self, names):
        # GH#17250 make sure result dtype is correct
        # GH#19043 make sure names are propogated correctly
        tdi = pd.TimedeltaIndex(['0 days', '1 day'], name=names[0])
        ser = Series([Timedelta(hours=3), Timedelta(hours=4)], name=names[1])
        expected = Series([Timedelta(hours=3), Timedelta(days=1, hours=4)],
                          name=names[2])

        result = tdi + ser
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'timedelta64[ns]'

        result = ser + tdi
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'timedelta64[ns]'

        expected = Series([Timedelta(hours=-3), Timedelta(days=1, hours=-4)],
                          name=names[2])

        result = tdi - ser
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'timedelta64[ns]'

        result = ser - tdi
        tm.assert_series_equal(result, -expected)
        assert result.dtype == 'timedelta64[ns]'

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_tdi_mul_int_series(self, names):
        # GH#19042
        tdi = pd.TimedeltaIndex(['0days', '1day', '2days', '3days', '4days'],
                                name=names[0])
        ser = Series([0, 1, 2, 3, 4], dtype=np.int64, name=names[1])

        expected = Series(['0days', '1day', '4days', '9days', '16days'],
                          dtype='timedelta64[ns]',
                          name=names[2])

        result = ser * tdi
        tm.assert_series_equal(result, expected)

        # The direct operation tdi * ser still needs to be fixed.
        result = ser.__rmul__(tdi)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_float_series_rdiv_tdi(self, names):
        # GH#19042
        # TODO: the direct operation TimedeltaIndex / Series still
        # needs to be fixed.
        tdi = pd.TimedeltaIndex(['0days', '1day', '2days', '3days', '4days'],
                                name=names[0])
        ser = Series([1.5, 3, 4.5, 6, 7.5], dtype=np.float64, name=names[1])

        expected = Series([tdi[n] / ser[n] for n in range(len(ser))],
                          dtype='timedelta64[ns]',
                          name=names[2])

        result = ser.__rdiv__(tdi)
        tm.assert_series_equal(result, expected)
