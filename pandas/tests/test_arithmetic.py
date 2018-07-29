# -*- coding: utf-8 -*-
# Arithmetc tests for DataFrame/Series/Index/Array classes that should
# behave identically.
from datetime import timedelta
import operator

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm

from pandas.core import ops
from pandas import (
    Timedelta, Timestamp, NaT, Series, TimedeltaIndex, DatetimeIndex)


def maybe_upcast_box(expected, right, box):
    # Wrap `expected` in the "higher" type between `box` and `type(right)`,
    # e.g. if box == pd.Index and right is a Series, we expected should
    # be wrapped in a Series
    if box is pd.DataFrame or isinstance(right, pd.DataFrame):
        return pd.DataFrame(expected)
    elif box is pd.Series or isinstance(right, pd.Series):
        return pd.Series(expected)
    else:
        # TODO: Once EA tests are in here, Index will be considered
        # "higher" than Array
        return pd.Index(expected)

# ------------------------------------------------------------------
# Fixtures

@pytest.fixture
def tdser():
    """
    Return a Series with dtype='timedelta64[ns]', including a NaT.
    """
    return Series(['59 Days', '59 Days', 'NaT'], dtype='timedelta64[ns]')


# ------------------------------------------------------------------
# Numeric dtypes Arithmetic with Timedelta Scalar

class TestNumericArraylikeArithmeticWithTimedeltaScalar(object):

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="block.eval incorrect",
                                             strict=True))
    ])
    @pytest.mark.parametrize('index', [
        pd.Int64Index(range(1, 11)),
        pd.UInt64Index(range(1, 11)),
        pd.Float64Index(range(1, 11)),
        pd.RangeIndex(1, 11)],
        ids=lambda x: type(x).__name__)
    @pytest.mark.parametrize('scalar_td', [
        Timedelta(days=1),
        Timedelta(days=1).to_timedelta64(),
        Timedelta(days=1).to_pytimedelta()],
        ids=lambda x: type(x).__name__)
    def test_index_mul_timedelta(self, scalar_td, index, box):
        # GH#19333

        if (box is Series and
                type(scalar_td) is timedelta and index.dtype == 'f8'):
            raise pytest.xfail(reason="Cannot multiply timedelta by float")

        expected = pd.timedelta_range('1 days', '10 days')

        index = tm.box_expected(index, box)
        expected = tm.box_expected(expected, box)

        result = index * scalar_td
        tm.assert_equal(result, expected)

        commute = scalar_td * index
        tm.assert_equal(commute, expected)

    @pytest.mark.parametrize('box', [pd.Index, Series, pd.DataFrame])
    @pytest.mark.parametrize('index', [
        pd.Int64Index(range(1, 3)),
        pd.UInt64Index(range(1, 3)),
        pd.Float64Index(range(1, 3)),
        pd.RangeIndex(1, 3)],
        ids=lambda x: type(x).__name__)
    @pytest.mark.parametrize('scalar_td', [
        Timedelta(days=1),
        Timedelta(days=1).to_timedelta64(),
        Timedelta(days=1).to_pytimedelta()],
        ids=lambda x: type(x).__name__)
    def test_index_rdiv_timedelta(self, scalar_td, index, box):

        if box is Series and type(scalar_td) is timedelta:
            raise pytest.xfail(reason="TODO: Figure out why this case fails")
        if box is pd.DataFrame and isinstance(scalar_td, timedelta):
            raise pytest.xfail(reason="TODO: Figure out why this case fails")

        expected = TimedeltaIndex(['1 Day', '12 Hours'])

        index = tm.box_expected(index, box)
        expected = tm.box_expected(expected, box)

        result = scalar_td / index
        tm.assert_equal(result, expected)

        with pytest.raises(TypeError):
            index / scalar_td


# ------------------------------------------------------------------
# Timedelta64[ns] dtype Arithmetic Operations

class TestTimedeltaArraylikeAddSubOps(object):
    # Tests for timedelta64[ns] __add__, __sub__, __radd__, __rsub__

    # -------------------------------------------------------------
    # Invalid Operations

    def test_tdi_add_str_invalid(self):
        # GH#13624
        tdi = TimedeltaIndex(['1 day', '2 days'])

        with pytest.raises(TypeError):
            tdi + 'a'
        with pytest.raises(TypeError):
            'a' + tdi

    @pytest.mark.parametrize('other', [3.14, np.array([2.0, 3.0])])
    @pytest.mark.parametrize('op', [operator.add, ops.radd,
                                    operator.sub, ops.rsub])
    def test_tdi_add_sub_float(self, op, other):
        tdi = TimedeltaIndex(['-1 days', '-1 days'])
        with pytest.raises(TypeError):
            op(tdi, other)

    @pytest.mark.parametrize('freq', [None, 'H'])
    def test_tdi_sub_period(self, freq):
        # GH#13078
        # not supported, check TypeError
        p = pd.Period('2011-01-01', freq='D')
        idx = TimedeltaIndex(['1 hours', '2 hours'], freq=freq)

        with pytest.raises(TypeError):
            idx - p

        with pytest.raises(TypeError):
            p - idx

    @pytest.mark.parametrize('op', [operator.add, ops.radd,
                                    operator.sub, ops.rsub])
    @pytest.mark.parametrize('pi_freq', ['D', 'W', 'Q', 'H'])
    @pytest.mark.parametrize('tdi_freq', [None, 'H'])
    def test_td_sub_pi(self, tdi_freq, pi_freq, op):
        # GH#20049 subtracting PeriodIndex should raise TypeError
        tdi = TimedeltaIndex(['1 hours', '2 hours'], freq=tdi_freq)
        dti = Timestamp('2018-03-07 17:16:40') + tdi
        pi = dti.to_period(pi_freq)
        with pytest.raises(TypeError):
            op(dti, pi)

    # -------------------------------------------------------------
    # Binary operations TimedeltaIndex and datetime-like

    def test_tdi_sub_timestamp_raises(self):
        idx = TimedeltaIndex(['1 day', '2 day'])
        msg = "cannot subtract a datelike from a TimedeltaIndex"
        with tm.assert_raises_regex(TypeError, msg):
            idx - Timestamp('2011-01-01')

    def test_tdi_add_timestamp(self):
        idx = TimedeltaIndex(['1 day', '2 day'])

        result = idx + Timestamp('2011-01-01')
        expected = DatetimeIndex(['2011-01-02', '2011-01-03'])
        tm.assert_index_equal(result, expected)

    def test_tdi_radd_timestamp(self):
        idx = TimedeltaIndex(['1 day', '2 day'])

        result = Timestamp('2011-01-01') + idx
        expected = DatetimeIndex(['2011-01-02', '2011-01-03'])
        tm.assert_index_equal(result, expected)

    # ------------------------------------------------------------------
    # Operations with int-like others

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

    def test_td64_series_add_intlike(self):
        # GH#19123
        tdi = TimedeltaIndex(['59 days', '59 days', 'NaT'])
        ser = Series(tdi)

        other = Series([20, 30, 40], dtype='uint8')

        with pytest.raises(TypeError):
            ser + 1
        with pytest.raises(TypeError):
            ser - 1

        with pytest.raises(TypeError):
            ser + other
        with pytest.raises(TypeError):
            ser - other

        with pytest.raises(TypeError):
            ser + other.values
        with pytest.raises(TypeError):
            ser - other.values

        with pytest.raises(TypeError):
            ser + pd.Index(other)
        with pytest.raises(TypeError):
            ser - pd.Index(other)

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
        Series([1, 2, 3])
    ])
    def test_td64series_add_sub_numeric_array_invalid(self, vector,
                                                      dtype, tdser):
        if type(vector) is Series and not dtype.startswith('float'):
            pytest.xfail(reason='GH#19123 integer interpreted as nanos')
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
    # Operations with datetime-like others

    def test_td64series_add_sub_timestamp(self):
        # GH#11925
        tdser = Series(pd.timedelta_range('1 day', periods=3))
        ts = Timestamp('2012-01-01')
        expected = Series(pd.date_range('2012-01-02', periods=3))
        tm.assert_series_equal(ts + tdser, expected)
        tm.assert_series_equal(tdser + ts, expected)

        expected2 = Series(pd.date_range('2011-12-31',
                                         periods=3, freq='-1D'))
        tm.assert_series_equal(ts - tdser, expected2)
        tm.assert_series_equal(ts + (-tdser), expected2)

        with pytest.raises(TypeError):
            tdser - ts

    # ------------------------------------------------------------------
    # Operations with timedelta-like others (including DateOffsets)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_td64_series_with_tdi(self, names):
        # GH#17250 make sure result dtype is correct
        # GH#19043 make sure names are propagated correctly
        tdi = TimedeltaIndex(['0 days', '1 day'], name=names[0])
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

    def test_td64_sub_NaT(self):
        # GH#18808
        ser = Series([NaT, Timedelta('1s')])
        res = ser - NaT
        expected = Series([NaT, NaT], dtype='timedelta64[ns]')
        tm.assert_series_equal(res, expected)


class TestTimedeltaArraylikeMulDivOps(object):
    # Tests for timedelta64[ns]
    # __mul__, __rmul__, __div__, __rdiv__, __floordiv__, __rfloordiv__

    # ------------------------------------------------------------------
    # __floordiv__, __rfloordiv__

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Incorrectly returns "
                                                    "m8[ns] instead of f8",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_td64arr_floordiv_tdscalar(self, box, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        expected = Series([0, 0, np.nan])

        td1 = tm.box_expected(td1, box)
        expected = tm.box_expected(expected, box)

        result = td1 // scalar_td
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Incorrectly casts to f8",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
    def test_td64arr_rfloordiv_tdscalar(self, box, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        expected = Series([1, 1, np.nan])

        td1 = tm.box_expected(td1, box)
        expected = tm.box_expected(expected, box)

        result = scalar_td // td1
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Returns m8[ns] dtype "
                                                    "instead of f8",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
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

    # ------------------------------------------------------------------
    # Operations with timedelta-like others

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="__mul__ op treats "
                                                    "timedelta other as i8; "
                                                    "rmul OK",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
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

    # ------------------------------------------------------------------
    # Operations with numeric others

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Returns object-dtype",
                                             strict=True))
    ], ids=lambda x: x.__name__)
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

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Returns object-dtype",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('two', [2, 2.0, np.array(2), np.array(2.0)])
    def test_td64arr_div_numeric_scalar(self, box, two, tdser):
        # GH#4521
        # divide/multiply by integers
        expected = Series(['29.5D', '29.5D', 'NaT'], dtype='timedelta64[ns]')

        tdser = tm.box_expected(tdser, box)
        expected = tm.box_expected(expected, box)

        result = tdser / two
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="broadcasts along "
                                                    "wrong axis",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [np.array([20, 30, 40]),
                                        pd.Index([20, 30, 40]),
                                        Series([20, 30, 40])])
    def test_td64arr_mul_numeric_array(self, box, vector, dtype, tdser):
        # GH#4521
        # divide/multiply by integers
        vector = vector.astype(dtype)

        expected = Series(['1180 Days', '1770 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        tdser = tm.box_expected(tdser, box)
        # TODO: Make this up-casting more systematic?
        box = Series if (box is pd.Index and type(vector) is Series) else box
        expected = tm.box_expected(expected, box)

        result = tdser * vector
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="broadcasts along "
                                                    "wrong axis",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [np.array([20, 30, 40]),
                                        pd.Index([20, 30, 40]),
                                        Series([20, 30, 40])],
                             ids=lambda x: type(x).__name__)
    def test_td64arr_rmul_numeric_array(self, box, vector, dtype, tdser):
        # GH#4521
        # divide/multiply by integers
        vector = vector.astype(dtype)

        expected = Series(['1180 Days', '1770 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        tdser = tm.box_expected(tdser, box)
        box = Series if (box is pd.Index and type(vector) is Series) else box
        expected = tm.box_expected(expected, box)

        result = vector * tdser
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="broadcasts along "
                                                    "wrong axis",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [np.array([20, 30, 40]),
                                        pd.Index([20, 30, 40]),
                                        Series([20, 30, 40])])
    def test_td64arr_div_numeric_array(self, box, vector, dtype, tdser):
        # GH#4521
        # divide/multiply by integers
        vector = vector.astype(dtype)
        expected = Series(['2.95D', '1D 23H 12m', 'NaT'],
                          dtype='timedelta64[ns]')

        tdser = tm.box_expected(tdser, box)
        box = Series if (box is pd.Index and type(vector) is Series) else box
        expected = tm.box_expected(expected, box)

        result = tdser / vector
        tm.assert_equal(result, expected)

        with pytest.raises(TypeError):
            vector / tdser

    # TODO: Should we be parametrizing over types for `ser` too?
    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="broadcasts along "
                                                    "wrong axis",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_td64arr_mul_int_series(self, box, names):
        # GH#19042 test for correct name attachment
        tdi = TimedeltaIndex(['0days', '1day', '2days', '3days', '4days'],
                             name=names[0])
        ser = Series([0, 1, 2, 3, 4], dtype=np.int64, name=names[1])

        expected = Series(['0days', '1day', '4days', '9days', '16days'],
                          dtype='timedelta64[ns]',
                          name=names[2])

        tdi = tm.box_expected(tdi, box)
        box = Series if (box is pd.Index and type(ser) is Series) else box
        expected = tm.box_expected(expected, box)

        result = ser * tdi
        tm.assert_equal(result, expected)

        # The direct operation tdi * ser still needs to be fixed.
        result = ser.__rmul__(tdi)
        tm.assert_equal(result, expected)

    # TODO: Should we be parametrizing over types for `ser` too?
    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pd.DataFrame
    ], ids=lambda x: x.__name__)
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

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="raises ValueError "
                                                    "instead of TypeError",
                                             strict=True))
    ])
    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=5, seconds=4),
        Timedelta('5m4s'),
        Timedelta('5m4s').to_timedelta64()])
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
