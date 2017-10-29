""" Test Matrix for arithmetic operations on DatetimeIndex, TimedeltaIndex,
and PeriodIndex
"""

import pytest

import pandas as pd
from pandas import Timestamp, Timedelta, NaT


class TestDatetimeLikeIndexArithmetic(object):
    # GH17991 checking for overflows and NaT masking on arithmetic ops

    # TODO: Fill out the matrix of allowed arithmetic operations:
    #   - __rsub__, __radd__
    #   - ops with scalars boxed in Index/Series/DataFrame/np.array
    #   - ops with scalars:
    #              NaT, Timestamp.min/max, Timedelta.min/max
    #              datetime, timedelta, date(?),
    #              relativedelta,
    #              np.datetime64, np.timedelta64,
    #              DateOffset,
    #              Period
    #   - timezone-aware variants
    #   - object-dtype, categorical dtype
    #   - PeriodIndex
    #   - consistency with .map(...) ?

    def test_timedeltaindex_add_timestamp_nat_masking(self):
        tdinat = pd.to_timedelta(['24658 days 11:15:00', 'NaT'])

        # tsneg.value < 0, tspos.value > 0
        tsneg = Timestamp('1950-01-01')
        tspos = Timestamp('1980-01-01')

        res1 = tdinat + tsneg
        assert res1[1] is NaT
        res2 = tdinat + tspos
        assert res2[1] is NaT

    def test_timedeltaindex_sub_timestamp_nat_masking(self):
        tdinat = pd.to_timedelta(['24658 days 11:15:00', 'NaT'])

        # tsneg.value < 0, tspos.value > 0
        tsneg = Timestamp('1950-01-01')
        tspos = Timestamp('1980-01-01')

        res1 = tdinat - tsneg
        assert res1[1] is NaT
        res2 = tdinat - tspos
        assert res2[1] is NaT

    def test_timedeltaindex_add_timestamp_overflow(self):
        tdimax = pd.to_timedelta(['24658 days 11:15:00', Timedelta.max])
        tdimin = pd.to_timedelta(['24658 days 11:15:00', Timedelta.min])

        # tsneg.value < 0, tspos.value > 0
        tsneg = Timestamp('1950-01-01')
        tspos = Timestamp('1980-01-01')

        res1 = tdimax + tsneg
        assert res1[1].value == Timedelta.max.value + tsneg.value
        res2 = tdimin + tspos
        assert res2[1].value == Timedelta.min.value + tspos.value

        with pytest.raises(OverflowError):
            tdimax + tspos

        with pytest.raises(OverflowError):
            tdimin + tsneg

    def test_timedeltaindex_add_timedelta_overflow(self):
        tdimax = pd.to_timedelta(['24658 days 11:15:00', Timedelta.max])
        tdimin = pd.to_timedelta(['24658 days 11:15:00', Timedelta.min])

        # tdpos.value > 0, tdneg.value < 0
        tdpos = Timedelta('1h')
        tdneg = Timedelta('-1h')

        with pytest.raises(OverflowError):
            tdimax + tdpos

        res2 = tdimax + tdneg
        assert res2[1].value == Timedelta.max.value + tdneg.value
        res3 = tdimin + tdpos
        assert res3[1].value == Timedelta.min.value + tdpos.value

        with pytest.raises(OverflowError):
            tdimin + tdneg

    def test_timedeltaindex_sub_timedelta_overflow(self):
        tdimax = pd.to_timedelta(['24658 days 11:15:00', Timedelta.max])
        tdimin = pd.to_timedelta(['24658 days 11:15:00', Timedelta.min])

        # tdpos.value > 0, tdneg.value < 0
        tdpos = Timedelta('1h')
        tdneg = Timedelta('-1h')

        res1 = tdimax - tdpos
        assert res1[1].value == Timedelta.max.value - tdpos.value

        with pytest.raises(OverflowError):
            tdimax - tdneg

        with pytest.raises(OverflowError):
            tdimin - tdpos

        res4 = tdimin - tdneg
        assert res4[1].value == Timedelta.min.value - tdneg.value

    def test_datetimeindex_add_nat_masking(self):
        # Checking for NaTs and checking that we don't get an OverflowError
        dtinat = pd.to_datetime(['now', 'NaT'])

        # tdpos.value > 0, tdneg.value < 0
        tdpos = Timedelta('1h')
        tdneg = Timedelta('-1h')

        res1 = dtinat + tdpos
        assert res1[1] is NaT
        res2 = dtinat + tdneg
        assert res2[1] is NaT

    def test_datetimeindex_sub_nat_masking(self):
        # Checking for NaTs and checking that we don't get an OverflowError
        dtinat = pd.to_datetime(['now', 'NaT'])

        # tdpos.value > 0, tdneg.value < 0
        tdpos = Timedelta('1h')
        tdneg = Timedelta('-1h')

        res1 = dtinat - tdpos
        assert res1[1] is NaT
        res2 = dtinat - tdneg
        assert res2[1] is NaT

    def test_datetimeindex_add_timedelta_overflow(self):
        dtimax = pd.to_datetime(['now', Timestamp.max])
        dtimin = pd.to_datetime(['now', Timestamp.min])

        # tdpos.value < 0, tdneg.value > 0
        tdpos = Timedelta('1h')
        tdneg = Timedelta('-1h')

        with pytest.raises(OverflowError):
            dtimax + tdpos

        res2 = dtimax + tdneg
        assert res2[1].value == Timestamp.max.value + tdneg.value

        res3 = dtimin + tdpos
        assert res3[1].value == Timestamp.min.value + tdpos.value

        with pytest.raises(OverflowError):
            dtimin + tdneg

    def test_datetimeindex_sub_timedelta_overflow(self):
        dtimax = pd.to_datetime(['now', Timestamp.max])
        dtimin = pd.to_datetime(['now', Timestamp.min])

        # tdpos.value < 0, tdneg.value > 0
        tdpos = Timedelta('1h')
        tdneg = Timedelta('-1h')

        res1 = dtimax - tdpos
        assert res1[1].value == Timestamp.max.value - tdpos.value

        with pytest.raises(OverflowError):
            dtimax - tdneg

        with pytest.raises(OverflowError):
            dtimin - tdpos

        res4 = dtimin - tdneg
        assert res4[1].value == Timestamp.min.value - tdneg.value

    def test_datetimeindex_sub_timestamp_overflow(self):
        dtimax = pd.to_datetime(['now', Timestamp.max])
        dtimin = pd.to_datetime(['now', Timestamp.min])

        # tsneg.value < 0, tspos.value > 0
        tsneg = Timestamp('1950-01-01')
        tspos = Timestamp('1980-01-01')

        with pytest.raises(OverflowError):
            dtimax - tsneg

        res2 = dtimax - tspos
        assert res2[1].value == Timestamp.max.value - tspos.value

        res3 = dtimin - tsneg
        assert res3[1].value == Timestamp.min.value - tsneg.value

        with pytest.raises(OverflowError):
            dtimin - tspos
