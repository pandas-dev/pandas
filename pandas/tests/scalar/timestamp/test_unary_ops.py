# -*- coding: utf-8 -*-
from datetime import datetime

import pytest
import pytz
from pytz import utc

import pandas.util.testing as tm
import pandas.util._test_decorators as td

from pandas.compat import PY3
from pandas._libs.tslibs.frequencies import _INVALID_FREQ_ERROR
from pandas import Timestamp


class TestTimestampUnaryOps(object):

    # --------------------------------------------------------------
    # Timestamp.round

    def test_round_day_naive(self):
        dt = Timestamp('20130101 09:10:11')
        result = dt.round('D')
        expected = Timestamp('20130101')
        assert result == expected

        dt = Timestamp('20130101 19:10:11')
        result = dt.round('D')
        expected = Timestamp('20130102')
        assert result == expected

        dt = Timestamp('20130201 12:00:00')
        result = dt.round('D')
        expected = Timestamp('20130202')
        assert result == expected

        dt = Timestamp('20130104 12:00:00')
        result = dt.round('D')
        expected = Timestamp('20130105')
        assert result == expected

    def test_round_tzaware(self):
        dt = Timestamp('20130101 09:10:11', tz='US/Eastern')
        result = dt.round('D')
        expected = Timestamp('20130101', tz='US/Eastern')
        assert result == expected

        dt = Timestamp('20130101 09:10:11', tz='US/Eastern')
        result = dt.round('s')
        assert result == dt

    def test_round_30min(self):
        # round
        dt = Timestamp('20130104 12:32:00')
        result = dt.round('30Min')
        expected = Timestamp('20130104 12:30:00')
        assert result == expected

    def test_round_subsecond(self):
        # GH#14440 & GH#15578
        result = Timestamp('2016-10-17 12:00:00.0015').round('ms')
        expected = Timestamp('2016-10-17 12:00:00.002000')
        assert result == expected

        result = Timestamp('2016-10-17 12:00:00.00149').round('ms')
        expected = Timestamp('2016-10-17 12:00:00.001000')
        assert result == expected

        ts = Timestamp('2016-10-17 12:00:00.0015')
        for freq in ['us', 'ns']:
            assert ts == ts.round(freq)

        result = Timestamp('2016-10-17 12:00:00.001501031').round('10ns')
        expected = Timestamp('2016-10-17 12:00:00.001501030')
        assert result == expected

    def test_round_nonstandard_freq(self):
        with tm.assert_produces_warning():
            Timestamp('2016-10-17 12:00:00.001501031').round('1010ns')

    def test_round_invalid_arg(self):
        stamp = Timestamp('2000-01-05 05:09:15.13')
        with tm.assert_raises_regex(ValueError, _INVALID_FREQ_ERROR):
            stamp.round('foo')

    @pytest.mark.parametrize('freq, expected', [
        ('D', Timestamp('2000-01-05 00:00:00')),
        ('H', Timestamp('2000-01-05 05:00:00')),
        ('S', Timestamp('2000-01-05 05:09:15'))])
    def test_round_frequencies(self, freq, expected):
        stamp = Timestamp('2000-01-05 05:09:15.13')

        result = stamp.round(freq=freq)
        assert result == expected

    def test_ceil(self):
        dt = Timestamp('20130101 09:10:11')
        result = dt.ceil('D')
        expected = Timestamp('20130102')
        assert result == expected

    def test_floor(self):
        dt = Timestamp('20130101 09:10:11')
        result = dt.floor('D')
        expected = Timestamp('20130101')
        assert result == expected

    # --------------------------------------------------------------
    # Timestamp.replace
    timezones = ['UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/US/Pacific']

    def test_replace_naive(self):
        # GH#14621, GH#7825
        ts = Timestamp('2016-01-01 09:00:00')
        result = ts.replace(hour=0)
        expected = Timestamp('2016-01-01 00:00:00')
        assert result == expected

    @pytest.mark.parametrize('tz', timezones)
    def test_replace_aware(self, tz):
        # GH#14621, GH#7825
        # replacing datetime components with and w/o presence of a timezone
        ts = Timestamp('2016-01-01 09:00:00', tz=tz)
        result = ts.replace(hour=0)
        expected = Timestamp('2016-01-01 00:00:00', tz=tz)
        assert result == expected

    @pytest.mark.parametrize('tz', timezones)
    def test_replace_preserves_nanos(self, tz):
        # GH#14621, GH#7825
        ts = Timestamp('2016-01-01 09:00:00.000000123', tz=tz)
        result = ts.replace(hour=0)
        expected = Timestamp('2016-01-01 00:00:00.000000123', tz=tz)
        assert result == expected

    @pytest.mark.parametrize('tz', timezones)
    def test_replace_multiple(self, tz):
        # GH#14621, GH#7825
        # replacing datetime components with and w/o presence of a timezone
        # test all
        ts = Timestamp('2016-01-01 09:00:00.000000123', tz=tz)
        result = ts.replace(year=2015, month=2, day=2, hour=0, minute=5,
                            second=5, microsecond=5, nanosecond=5)
        expected = Timestamp('2015-02-02 00:05:05.000005005', tz=tz)
        assert result == expected

    @pytest.mark.parametrize('tz', timezones)
    def test_replace_invalid_kwarg(self, tz):
        # GH#14621, GH#7825
        ts = Timestamp('2016-01-01 09:00:00.000000123', tz=tz)
        with pytest.raises(TypeError):
            ts.replace(foo=5)

    @pytest.mark.parametrize('tz', timezones)
    def test_replace_integer_args(self, tz):
        # GH#14621, GH#7825
        ts = Timestamp('2016-01-01 09:00:00.000000123', tz=tz)
        with pytest.raises(ValueError):
            ts.replace(hour=0.1)

    def test_replace_tzinfo_equiv_tz_localize_none(self):
        # GH#14621, GH#7825
        # assert conversion to naive is the same as replacing tzinfo with None
        ts = Timestamp('2013-11-03 01:59:59.999999-0400', tz='US/Eastern')
        assert ts.tz_localize(None) == ts.replace(tzinfo=None)

    @td.skip_if_windows
    def test_replace_tzinfo(self):
        # GH#15683
        dt = datetime(2016, 3, 27, 1)
        tzinfo = pytz.timezone('CET').localize(dt, is_dst=False).tzinfo

        result_dt = dt.replace(tzinfo=tzinfo)
        result_pd = Timestamp(dt).replace(tzinfo=tzinfo)

        if PY3:
            # datetime.timestamp() converts in the local timezone
            with tm.set_timezone('UTC'):
                assert result_dt.timestamp() == result_pd.timestamp()

        assert result_dt == result_pd
        assert result_dt == result_pd.to_pydatetime()

        result_dt = dt.replace(tzinfo=tzinfo).replace(tzinfo=None)
        result_pd = Timestamp(dt).replace(tzinfo=tzinfo).replace(tzinfo=None)

        if PY3:
            # datetime.timestamp() converts in the local timezone
            with tm.set_timezone('UTC'):
                assert result_dt.timestamp() == result_pd.timestamp()

        assert result_dt == result_pd
        assert result_dt == result_pd.to_pydatetime()

    # --------------------------------------------------------------

    @td.skip_if_windows
    def test_timestamp(self):
        # GH#17329
        # tz-naive --> treat it as if it were UTC for purposes of timestamp()
        ts = Timestamp.now()
        uts = ts.replace(tzinfo=utc)
        assert ts.timestamp() == uts.timestamp()

        tsc = Timestamp('2014-10-11 11:00:01.12345678', tz='US/Central')
        utsc = tsc.tz_convert('UTC')

        # utsc is a different representation of the same time
        assert tsc.timestamp() == utsc.timestamp()

        if PY3:
            # datetime.timestamp() converts in the local timezone
            with tm.set_timezone('UTC'):

                # should agree with datetime.timestamp method
                dt = ts.to_pydatetime()
                assert dt.timestamp() == ts.timestamp()
