# -*- coding: utf-8 -*-
"""
Tests for Timestamp timezone-related methods
"""

import pytest
from pytz.exceptions import AmbiguousTimeError, NonExistentTimeError

import pandas.util.testing as tm
from pandas import Timestamp, NaT


class TestTimestampTZOperations(object):
    # --------------------------------------------------------------
    # Timestamp.tz_localize

    def test_tz_localize_ambiguous(self):
        ts = Timestamp('2014-11-02 01:00')
        ts_dst = ts.tz_localize('US/Eastern', ambiguous=True)
        ts_no_dst = ts.tz_localize('US/Eastern', ambiguous=False)

        assert (ts_no_dst.value - ts_dst.value) / 1e9 == 3600
        with pytest.raises(ValueError):
            ts.tz_localize('US/Eastern', ambiguous='infer')

        # GH#8025
        with tm.assert_raises_regex(TypeError,
                                    'Cannot localize tz-aware Timestamp, '
                                    'use tz_convert for conversions'):
            Timestamp('2011-01-01', tz='US/Eastern').tz_localize('Asia/Tokyo')

        with tm.assert_raises_regex(TypeError,
                                    'Cannot convert tz-naive Timestamp, '
                                    'use tz_localize to localize'):
            Timestamp('2011-01-01').tz_convert('Asia/Tokyo')

    @pytest.mark.parametrize('stamp, tz', [
        ('2015-03-08 02:00', 'US/Eastern'),
        ('2015-03-08 02:30', 'US/Pacific'),
        ('2015-03-29 02:00', 'Europe/Paris'),
        ('2015-03-29 02:30', 'Europe/Belgrade')])
    def test_tz_localize_nonexistent(self, stamp, tz):
        # GH#13057
        ts = Timestamp(stamp)
        with pytest.raises(NonExistentTimeError):
            ts.tz_localize(tz)
        with pytest.raises(NonExistentTimeError):
            ts.tz_localize(tz, errors='raise')
        assert ts.tz_localize(tz, errors='coerce') is NaT

    def test_tz_localize_errors_ambiguous(self):
        # GH#13057
        ts = Timestamp('2015-11-1 01:00')
        with pytest.raises(AmbiguousTimeError):
            ts.tz_localize('US/Pacific', errors='coerce')

    @pytest.mark.parametrize('tz', ['UTC', 'Asia/Tokyo',
                                    'US/Eastern', 'dateutil/US/Pacific'])
    @pytest.mark.parametrize('stamp', ['2014-02-01 09:00', '2014-07-08 09:00',
                                       '2014-11-01 17:00', '2014-11-05 00:00'])
    def test_tz_localize_roundtrip(self, stamp, tz):
        ts = Timestamp(stamp)
        localized = ts.tz_localize(tz)
        assert localized == Timestamp(stamp, tz=tz)

        with pytest.raises(TypeError):
            localized.tz_localize(tz)

        reset = localized.tz_localize(None)
        assert reset == ts
        assert reset.tzinfo is None

    # ------------------------------------------------------------------
    # Timestamp.tz_convert

    @pytest.mark.parametrize('tz', ['UTC', 'Asia/Tokyo',
                                    'US/Eastern', 'dateutil/US/Pacific'])
    @pytest.mark.parametrize('stamp', ['2014-02-01 09:00', '2014-07-08 09:00',
                                       '2014-11-01 17:00', '2014-11-05 00:00'])
    def test_tz_convert_roundtrip(self, stamp, tz):
        ts = Timestamp(stamp, tz='UTC')
        converted = ts.tz_convert(tz)

        reset = converted.tz_convert(None)
        assert reset == Timestamp(stamp)
        assert reset.tzinfo is None
        assert reset == converted.tz_convert('UTC').tz_localize(None)
