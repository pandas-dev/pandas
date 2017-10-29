""" Test Matrix for arithmetic operations on DatetimeIndex, TimedeltaIndex,
and PeriodIndex
"""
from datetime import datetime, timedelta

import pytest
import numpy as np

import pandas as pd
from pandas import Timestamp, Timedelta, NaT

tdinat = pd.to_timedelta(['24658 days 11:15:00', 'NaT'])
tdimax = pd.to_timedelta(['24658 days 11:15:00', Timedelta.max])
tdimin = pd.to_timedelta(['24658 days 11:15:00', Timedelta.min])

dtinat = pd.to_datetime(['now', 'NaT'])
dtimax = pd.to_datetime(['now', Timestamp.max])
dtimin = pd.to_datetime(['now', Timestamp.min])


tspos = Timestamp('1980-01-01')
ts_pos_variants = [tspos,
                   datetime(1980, 1, 1),
                   np.datetime64('1980-01-01').astype('M8[ns]'),
                   np.datetime64('1980-01-01').astype('M8[D]')]

tsneg = Timestamp('1950-01-01')
ts_neg_variants = [tsneg,
                   datetime(1950, 1, 1),
                   np.datetime64('1950-01-01').astype('M8[ns]'),
                   np.datetime64('1950-01-01').astype('M8[D]')]

tdpos = Timedelta('1h')
td_pos_variants = [tdpos,
                   tdpos.to_pytimedelta(),
                   tdpos.to_timedelta64()]

tdneg = Timedelta('-1h')
td_neg_variants = [tdneg,
                   tdneg.to_pytimedelta(),
                   tdneg.to_timedelta64()]


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
    #   - versions with near-min/max values

    def test_timedeltaindex_add_timestamp_nat_masking(self):
        for variant in ts_neg_variants:
            res = tdinat + variant
            assert res[1] is NaT

        for variant in ts_pos_variants:
            res = tdinat + variant
            assert res[1] is NaT

    def test_timedeltaindex_add_timestamp_overflow(self):
        expected = Timedelta.max.value + tsneg.value
        for variant in ts_neg_variants:
            res = tdimax + variant
            assert res[1].value == expected

        expected = Timedelta.min.value + tspos.value
        for variant in ts_pos_variants:
            res = tdimin + variant
            assert res[1].value == expected

        for variant in ts_pos_variants:
            with pytest.raises(OverflowError):
                tdimax + variant

        for variant in ts_neg_variants:
            with pytest.raises(OverflowError):
                tdimin + variant

    def test_timedeltaindex_add_timedelta_overflow(self):
        for variant in td_pos_variants:
            with pytest.raises(OverflowError):
                tdimax + variant

        expected = Timedelta.max.value + tdneg.value
        for variant in td_neg_variants:
            res = tdimax + variant
            assert res[1].value == expected

        expected = Timedelta.min.value + tdpos.value
        for variant in td_pos_variants:
            res = tdimin + variant
            assert res[1].value == expected

        for variant in td_neg_variants:
            with pytest.raises(OverflowError):
                tdimin + variant

    def test_timedeltaindex_sub_timedelta_overflow(self):
        expected = Timedelta.max.value - tdpos.value
        for variant in td_pos_variants:
            res1 = tdimax - variant
            assert res1[1].value == expected

        for variant in td_neg_variants:
            with pytest.raises(OverflowError):
                tdimax - variant

        for variant in td_pos_variants:
            with pytest.raises(OverflowError):
                tdimin - variant

        expected = Timedelta.min.value - tdneg.value
        for variant in td_neg_variants:
            res = tdimin - variant
            assert res[1].value == expected

    def test_datetimeindex_add_nat_masking(self):
        # Checking for NaTs and checking that we don't get an OverflowError
        for variant in td_pos_variants:
            res = dtinat + variant
            assert res[1] is NaT

        for variant in td_neg_variants:
            res = dtinat + variant
            assert res[1] is NaT

    def test_datetimeindex_sub_nat_masking(self):
        # Checking for NaTs and checking that we don't get an OverflowError
        for variant in td_pos_variants:
            res = dtinat - variant
            assert res[1] is NaT

        for variant in td_neg_variants:
            res = dtinat - variant
            assert res[1] is NaT

    def test_datetimeindex_add_timedelta_overflow(self):
        for variant in td_pos_variants:
            with pytest.raises(OverflowError):
                dtimax + variant

        expected = Timestamp.max.value + tdneg.value
        for variant in td_neg_variants:
            res = dtimax + variant
            assert res[1].value == expected

        expected = Timestamp.min.value + tdpos.value
        for variant in td_pos_variants:
            res = dtimin + variant
            assert res[1].value == expected

        for variant in td_neg_variants:
            with pytest.raises(OverflowError):
                dtimin + variant

    def test_datetimeindex_sub_timedelta_overflow(self):
        expected = Timestamp.max.value - tdpos.value
        for variant in td_pos_variants:
            res = dtimax - variant
            assert res[1].value == expected

        for variant in td_neg_variants:
            with pytest.raises(OverflowError):
                dtimax - variant

        for variant in td_pos_variants:
            with pytest.raises(OverflowError):
                dtimin - variant

        expected = Timestamp.min.value - tdneg.value
        for variant in td_neg_variants:
            res = dtimin - variant
            assert res[1].value == expected

    def test_datetimeindex_sub_timestamp_overflow(self):
        for variant in ts_neg_variants:
            with pytest.raises(OverflowError):
                dtimax - variant

        expected = Timestamp.max.value - tspos.value
        for variant in ts_pos_variants:
            res = dtimax - variant
            assert res[1].value == expected

        expected = Timestamp.min.value - tsneg.value
        for variant in ts_neg_variants:
            res = dtimin - variant
            assert res[1].value == expected

        for variant in ts_pos_variants:
            with pytest.raises(OverflowError):
                dtimin - variant
