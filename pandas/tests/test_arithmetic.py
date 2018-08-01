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
from pandas.errors import NullFrequencyError
from pandas._libs.tslibs import IncompatibleFrequency
from pandas import (
    timedelta_range,
    Timedelta, Timestamp, NaT, Series, TimedeltaIndex, DatetimeIndex)


# ------------------------------------------------------------------
# Fixtures

@pytest.fixture
def tdser():
    """
    Return a Series with dtype='timedelta64[ns]', including a NaT.
    """
    return Series(['59 Days', '59 Days', 'NaT'], dtype='timedelta64[ns]')


@pytest.fixture(params=[pd.offsets.Hour(2), timedelta(hours=2),
                        np.timedelta64(2, 'h'), Timedelta(hours=2)],
                ids=lambda x: type(x).__name__)
def delta(request):
    """
    Several ways of representing two hours
    """
    return request.param


@pytest.fixture(params=[timedelta(minutes=5, seconds=4),
                        Timedelta('5m4s'),
                        Timedelta('5m4s').to_timedelta64()],
                ids=lambda x: type(x).__name__)
def scalar_td(request):
    """
    Several variants of Timedelta scalars representing 5 minutes and 4 seconds
    """
    return request.param


@pytest.fixture(params=[pd.Index, Series, pd.DataFrame],
                ids=lambda x: x.__name__)
def box(request):
    """
    Several array-like containers that should have effectively identical
    behavior with respect to arithmetic operations.
    """
    return request.param


@pytest.fixture(params=[pd.Index,
                        Series,
                        pytest.param(pd.DataFrame,
                                     marks=pytest.mark.xfail(strict=True))],
                ids=lambda x: x.__name__)
def box_df_fail(request):
    """
    Fixture equivalent to `box` fixture but xfailing the DataFrame case.
    """
    return request.param


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
    def test_numeric_arr_mul_tdscalar(self, scalar_td, index, box):
        # GH#19333

        if (box is Series and
                type(scalar_td) is timedelta and index.dtype == 'f8'):
            raise pytest.xfail(reason="Cannot multiply timedelta by float")

        expected = timedelta_range('1 days', '10 days')

        index = tm.box_expected(index, box)
        expected = tm.box_expected(expected, box)

        result = index * scalar_td
        tm.assert_equal(result, expected)

        commute = scalar_td * index
        tm.assert_equal(commute, expected)

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
    def test_numeric_arr_rdiv_tdscalar(self, scalar_td, index, box):

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
        tdi = tm.box_expected(tdi, box)

        if box is pd.DataFrame and op in [operator.add, operator.sub]:
            pytest.xfail(reason="Tries to align incorrectly, "
                                "raises ValueError")

        with pytest.raises(TypeError):
            op(tdi, other)

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Tries to cast df to "
                                                    "Period",
                                             strict=True,
                                             raises=IncompatibleFrequency))
    ], ids=lambda x: x.__name__)
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

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="broadcasts along "
                                                    "wrong axis",
                                             raises=ValueError,
                                             strict=True))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('pi_freq', ['D', 'W', 'Q', 'H'])
    @pytest.mark.parametrize('tdi_freq', [None, 'H'])
    def test_td64arr_sub_pi(self, box, tdi_freq, pi_freq):
        # GH#20049 subtracting PeriodIndex should raise TypeError
        tdi = TimedeltaIndex(['1 hours', '2 hours'], freq=tdi_freq)
        dti = Timestamp('2018-03-07 17:16:40') + tdi
        pi = dti.to_period(pi_freq)

        # TODO: parametrize over box for pi?
        tdi = tm.box_expected(tdi, box)
        with pytest.raises(TypeError):
            tdi - pi

    # -------------------------------------------------------------
    # Binary operations td64 arraylike and datetime-like

    def test_td64arr_sub_timestamp_raises(self, box):
        idx = TimedeltaIndex(['1 day', '2 day'])
        idx = tm.box_expected(idx, box)

        msg = "cannot subtract a datelike from|Could not operate"
        with tm.assert_raises_regex(TypeError, msg):
            idx - Timestamp('2011-01-01')

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Returns object dtype",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    def test_td64arr_add_timestamp(self, box):
        idx = TimedeltaIndex(['1 day', '2 day'])
        expected = DatetimeIndex(['2011-01-02', '2011-01-03'])

        idx = tm.box_expected(idx, box)
        expected = tm.box_expected(expected, box)

        result = idx + Timestamp('2011-01-01')
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Returns object dtype",
                                             strict=True))
    ], ids=lambda x: x.__name__)
    def test_td64_radd_timestamp(self, box):
        idx = TimedeltaIndex(['1 day', '2 day'])
        expected = DatetimeIndex(['2011-01-02', '2011-01-03'])

        idx = tm.box_expected(idx, box)
        expected = tm.box_expected(expected, box)

        # TODO: parametrize over scalar datetime types?
        result = Timestamp('2011-01-01') + idx
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Returns object dtype "
                                                    "instead of "
                                                    "datetime64[ns]",
                                             strict=True))
    ], ids=lambda x: x.__name__)
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

    # ------------------------------------------------------------------
    # Operations with int-like others

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Attempts to broadcast "
                                                    "incorrectly",
                                             strict=True, raises=ValueError))
    ], ids=lambda x: x.__name__)
    def test_td64arr_add_int_series_invalid(self, box, tdser):
        tdser = tm.box_expected(tdser, box)
        err = TypeError if box is not pd.Index else NullFrequencyError
        with pytest.raises(err):
            tdser + Series([2, 3, 4])

    @pytest.mark.parametrize('box', [
        pd.Index,
        pytest.param(Series,
                     marks=pytest.mark.xfail(reason="GH#19123 integer "
                                                    "interpreted as "
                                                    "nanoseconds",
                                             strict=True)),
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Attempts to broadcast "
                                                    "incorrectly",
                                             strict=True, raises=ValueError))
    ], ids=lambda x: x.__name__)
    def test_td64arr_radd_int_series_invalid(self, box, tdser):
        tdser = tm.box_expected(tdser, box)
        err = TypeError if box is not pd.Index else NullFrequencyError
        with pytest.raises(err):
            Series([2, 3, 4]) + tdser

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Attempts to broadcast "
                                                    "incorrectly",
                                             strict=True, raises=ValueError))
    ], ids=lambda x: x.__name__)
    def test_td64arr_sub_int_series_invalid(self, box, tdser):
        tdser = tm.box_expected(tdser, box)
        err = TypeError if box is not pd.Index else NullFrequencyError
        with pytest.raises(err):
            tdser - Series([2, 3, 4])

    @pytest.mark.xfail(reason='GH#19123 integer interpreted as nanoseconds',
                       strict=True)
    def test_td64arr_rsub_int_series_invalid(self, box, tdser):
        tdser = tm.box_expected(tdser, box)
        with pytest.raises(TypeError):
            Series([2, 3, 4]) - tdser

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Tries to broadcast "
                                                    "incorrectly",
                                             strict=True, raises=ValueError))
    ], ids=lambda x: x.__name__)
    def test_td64arr_add_intlike(self, box):
        # GH#19123
        tdi = TimedeltaIndex(['59 days', '59 days', 'NaT'])
        ser = tm.box_expected(tdi, box)
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

        if box is pd.DataFrame and isinstance(scalar, np.ndarray):
            # raises ValueError
            pytest.xfail(reason="DataFrame to broadcast incorrectly")

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

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Tries to broadcast "
                                                    "incorrectly",
                                             strict=True, raises=ValueError))
    ], ids=lambda x: x.__name__)
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
        if type(vec) is Series and not dtype.startswith('float'):
            pytest.xfail(reason='GH#19123 integer interpreted as nanos')

        tdser = tm.box_expected(tdser, box)
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
    # Operations with timedelta-like others (including DateOffsets)

    # TODO: parametrize over [add, sub, radd, rsub]?
    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="Tries to broadcast "
                                                    "incorrectly leading "
                                                    "to alignment error",
                                             strict=True, raises=ValueError))
    ], ids=lambda x: x.__name__)
    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_td64arr_add_sub_tdi(self, box, names):
        # GH#17250 make sure result dtype is correct
        # GH#19043 make sure names are propagated correctly
        tdi = TimedeltaIndex(['0 days', '1 day'], name=names[0])
        ser = Series([Timedelta(hours=3), Timedelta(hours=4)], name=names[1])
        expected = Series([Timedelta(hours=3), Timedelta(days=1, hours=4)],
                          name=names[2])

        ser = tm.box_expected(ser, box)
        expected = tm.box_expected(expected, box)

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
        expected = tm.box_expected(expected, box)

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

    def test_td64arr_add_timedeltalike(self, delta, box):
        # only test adding/sub offsets as + is now numeric
        if box is pd.DataFrame and isinstance(delta, pd.DateOffset):
            pytest.xfail(reason="Returns object dtype instead of m8[ns]")

        rng = timedelta_range('1 days', '10 days')
        expected = timedelta_range('1 days 02:00:00', '10 days 02:00:00',
                                   freq='D')
        rng = tm.box_expected(rng, box)
        expected = tm.box_expected(expected, box)

        result = rng + delta
        tm.assert_equal(result, expected)

    def test_td64arr_sub_timedeltalike(self, delta, box):
        # only test adding/sub offsets as - is now numeric
        if box is pd.DataFrame and isinstance(delta, pd.DateOffset):
            pytest.xfail(reason="Returns object dtype instead of m8[ns]")

        rng = timedelta_range('1 days', '10 days')
        expected = timedelta_range('0 days 22:00:00', '9 days 22:00:00')

        rng = tm.box_expected(rng, box)
        expected = tm.box_expected(expected, box)

        result = rng - delta
        tm.assert_equal(result, expected)


class TestTimedeltaArraylikeMulDivOps(object):
    # Tests for timedelta64[ns]
    # __mul__, __rmul__, __div__, __rdiv__, __floordiv__, __rfloordiv__

    # ------------------------------------------------------------------
    # Multiplication
    # organized with scalar others first, then array-like

    def test_td64arr_mul_int(self, box_df_fail):
        box = box_df_fail  # DataFrame op returns object instead of m8[ns]
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        idx = tm.box_expected(idx, box)

        result = idx * 1
        tm.assert_equal(result, idx)

        result = 1 * idx
        tm.assert_equal(result, idx)

    def test_td64arr_mul_tdlike_scalar_raises(self, delta, box):
        if box is pd.DataFrame and not isinstance(delta, pd.DateOffset):
            pytest.xfail(reason="returns m8[ns] instead of raising")

        rng = timedelta_range('1 days', '10 days', name='foo')
        rng = tm.box_expected(rng, box)
        with pytest.raises(TypeError):
            rng * delta

    def test_tdi_mul_int_array_zerodim(self, box_df_fail):
        box = box_df_fail  # DataFrame op returns object dtype
        rng5 = np.arange(5, dtype='int64')
        idx = TimedeltaIndex(rng5)
        expected = TimedeltaIndex(rng5 * 5)

        idx = tm.box_expected(idx, box)
        expected = tm.box_expected(expected, box)

        result = idx * np.array(5, dtype='int64')
        tm.assert_equal(result, expected)

    def test_tdi_mul_int_array(self, box_df_fail):
        box = box_df_fail  # DataFrame tries to broadcast incorrectly
        rng5 = np.arange(5, dtype='int64')
        idx = TimedeltaIndex(rng5)
        expected = TimedeltaIndex(rng5 ** 2)

        idx = tm.box_expected(idx, box)
        expected = tm.box_expected(expected, box)

        result = idx * rng5
        tm.assert_equal(result, expected)

    def test_tdi_mul_int_series(self, box_df_fail):
        box = box_df_fail  # DataFrame tries to broadcast incorrectly
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        expected = TimedeltaIndex(np.arange(5, dtype='int64') ** 2)

        idx = tm.box_expected(idx, box)

        box2 = pd.Series if box is pd.Index else box
        expected = tm.box_expected(expected, box2)

        result = idx * pd.Series(np.arange(5, dtype='int64'))
        tm.assert_equal(result, expected)

    def test_tdi_mul_float_series(self, box_df_fail):
        box = box_df_fail  # DataFrame tries to broadcast incorrectly
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        idx = tm.box_expected(idx, box)

        rng5f = np.arange(5, dtype='float64')
        expected = TimedeltaIndex(rng5f * (rng5f + 0.1))
        box2 = pd.Series if box is pd.Index else box
        expected = tm.box_expected(expected, box2)

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
    def test_tdi_rmul_arraylike(self, other, box_df_fail):
        # RangeIndex fails to return NotImplemented, for others
        # DataFrame tries to broadcast incorrectly
        box = box_df_fail

        tdi = TimedeltaIndex(['1 Day'] * 10)
        expected = timedelta_range('1 days', '10 days')

        tdi = tm.box_expected(tdi, box)
        expected = tm.box_expected(expected, box)

        result = other * tdi
        tm.assert_equal(result, expected)
        commute = tdi * other
        tm.assert_equal(commute, expected)

    # ------------------------------------------------------------------
    # __div__

    def test_td64arr_div_nat_invalid(self, box_df_fail):
        # don't allow division by NaT (maybe could in the future)
        box = box_df_fail  # DataFrame returns all-NaT instead of raising
        rng = timedelta_range('1 days', '10 days', name='foo')
        rng = tm.box_expected(rng, box)
        with pytest.raises(TypeError):
            rng / pd.NaT

    def test_td64arr_div_int(self, box_df_fail):
        box = box_df_fail  # DataFrame returns object dtype instead of m8[ns]
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        idx = tm.box_expected(idx, box)

        result = idx / 1
        tm.assert_equal(result, idx)

    def test_tdi_div_tdlike_scalar(self, delta, box_df_fail):
        box = box_df_fail  # DataFrame op returns m8[ns] instead of float64
        rng = timedelta_range('1 days', '10 days', name='foo')
        expected = pd.Float64Index((np.arange(10) + 1) * 12, name='foo')

        rng = tm.box_expected(rng, box)
        expected = tm.box_expected(expected, box)

        result = rng / delta
        tm.assert_equal(result, expected)

    def test_tdi_div_tdlike_scalar_with_nat(self, delta, box_df_fail):
        box = box_df_fail  # DataFrame op returns m8[ns] instead of float64
        rng = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')
        expected = pd.Float64Index([12, np.nan, 24], name='foo')

        rng = tm.box_expected(rng, box)
        expected = tm.box_expected(expected, box)

        result = rng / delta
        tm.assert_equal(result, expected)

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

    def test_td64arr_floordiv_int(self, box_df_fail):
        box = box_df_fail  # DataFrame returns object dtype
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        idx = tm.box_expected(idx, box)
        result = idx // 1
        tm.assert_equal(result, idx)

    def test_td64arr_floordiv_tdlike_scalar(self, delta, box_df_fail):
        box = box_df_fail  # DataFrame returns m8[ns] instead of int64 dtype
        tdi = timedelta_range('1 days', '10 days', name='foo')
        expected = pd.Int64Index((np.arange(10) + 1) * 12, name='foo')

        tdi = tm.box_expected(tdi, box)
        expected = tm.box_expected(expected, box)

        result = tdi // delta
        tm.assert_equal(result, expected)

    # TODO: Is this redundant with test_td64arr_floordiv_tdlike_scalar?
    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=10, seconds=7),
        Timedelta('10m7s'),
        Timedelta('10m7s').to_timedelta64()
    ], ids=lambda x: type(x).__name__)
    def test_td64arr_rfloordiv_tdlike_scalar(self, scalar_td, box_df_fail):
        # GH#19125
        box = box_df_fail  # DataFrame op returns m8[ns] instead of f8 dtype
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

    @pytest.mark.parametrize('box', [
        pd.Index,
        Series,
        pytest.param(pd.DataFrame,
                     marks=pytest.mark.xfail(reason="__mul__ op treats "
                                                    "timedelta other as i8; "
                                                    "rmul OK",
                                             strict=True))
    ], ids=lambda x: x.__name__)
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
                                        Series([20, 30, 40])],
                             ids=lambda x: type(x).__name__)
    @pytest.mark.parametrize('op', [operator.mul, ops.rmul])
    def test_td64arr_rmul_numeric_array(self, op, box, vector, dtype, tdser):
        # GH#4521
        # divide/multiply by integers
        vector = vector.astype(dtype)

        expected = Series(['1180 Days', '1770 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        tdser = tm.box_expected(tdser, box)
        # TODO: Make this up-casting more systematic?
        box = Series if (box is pd.Index and type(vector) is Series) else box
        expected = tm.box_expected(expected, box)

        result = op(vector, tdser)
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
