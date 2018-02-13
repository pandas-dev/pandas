# -*- coding: utf-8 -*-

import pytest

from pandas.errors import OutOfBoundsDatetime
from pandas._libs.tslibs.frequencies import get_freq
from pandas._libs.tslibs.period import period_ordinal, period_asfreq, Period


class TestPeriodFreqConversion(object):
    @pytest.mark.parametrize('freq', ['A', 'Q', 'M', 'W', 'B', 'D'])
    def test_asfreq_near_zero(self, freq):
        # GH#19643, GH#19650
        per = Period('0001-01-01', freq=freq)
        tup1 = (per.year, per.hour, per.day)

        prev = per - 1
        assert (per - 1).ordinal == per.ordinal - 1
        tup2 = (prev.year, prev.month, prev.day)
        assert tup2 < tup1

    @pytest.mark.xfail(reason='GH#19643 period_helper asfreq functions fail '
                              'to check for overflows')
    def test_to_timestamp_out_of_bounds(self):
        # GH#19643, currently gives Timestamp('1754-08-30 22:43:41.128654848')
        per = Period('0001-01-01', freq='B')
        with pytest.raises(OutOfBoundsDatetime):
            per.to_timestamp()

    def test_intraday_conversion_factors(self):
        assert period_asfreq(1, get_freq('D'), get_freq('H'), False) == 24
        assert period_asfreq(1, get_freq('D'), get_freq('T'), False) == 1440
        assert period_asfreq(1, get_freq('D'), get_freq('S'), False) == 86400
        assert period_asfreq(1, get_freq('D'),
                             get_freq('L'), False) == 86400000
        assert period_asfreq(1, get_freq('D'),
                             get_freq('U'), False) == 86400000000
        assert period_asfreq(1, get_freq('D'),
                             get_freq('N'), False) == 86400000000000

        assert period_asfreq(1, get_freq('H'), get_freq('T'), False) == 60
        assert period_asfreq(1, get_freq('H'), get_freq('S'), False) == 3600
        assert period_asfreq(1, get_freq('H'),
                             get_freq('L'), False) == 3600000
        assert period_asfreq(1, get_freq('H'),
                             get_freq('U'), False) == 3600000000
        assert period_asfreq(1, get_freq('H'),
                             get_freq('N'), False) == 3600000000000

        assert period_asfreq(1, get_freq('T'), get_freq('S'), False) == 60
        assert period_asfreq(1, get_freq('T'), get_freq('L'), False) == 60000
        assert period_asfreq(1, get_freq('T'),
                             get_freq('U'), False) == 60000000
        assert period_asfreq(1, get_freq('T'),
                             get_freq('N'), False) == 60000000000

        assert period_asfreq(1, get_freq('S'), get_freq('L'), False) == 1000
        assert period_asfreq(1, get_freq('S'),
                             get_freq('U'), False) == 1000000
        assert period_asfreq(1, get_freq('S'),
                             get_freq('N'), False) == 1000000000

        assert period_asfreq(1, get_freq('L'), get_freq('U'), False) == 1000
        assert period_asfreq(1, get_freq('L'),
                             get_freq('N'), False) == 1000000

        assert period_asfreq(1, get_freq('U'), get_freq('N'), False) == 1000

    def test_period_ordinal_start_values(self):
        # information for 1.1.1970
        assert period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0, get_freq('A')) == 0
        assert period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0, get_freq('M')) == 0
        assert period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0, get_freq('W')) == 1
        assert period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0, get_freq('D')) == 0
        assert period_ordinal(1970, 1, 1, 0, 0, 0, 0, 0, get_freq('B')) == 0

    def test_period_ordinal_week(self):
        assert period_ordinal(1970, 1, 4, 0, 0, 0, 0, 0, get_freq('W')) == 1
        assert period_ordinal(1970, 1, 5, 0, 0, 0, 0, 0, get_freq('W')) == 2
        assert period_ordinal(2013, 10, 6, 0,
                              0, 0, 0, 0, get_freq('W')) == 2284
        assert period_ordinal(2013, 10, 7, 0,
                              0, 0, 0, 0, get_freq('W')) == 2285

    def test_period_ordinal_business_day(self):
        # Thursday
        assert period_ordinal(2013, 10, 3, 0,
                              0, 0, 0, 0, get_freq('B')) == 11415
        # Friday
        assert period_ordinal(2013, 10, 4, 0,
                              0, 0, 0, 0, get_freq('B')) == 11416
        # Saturday
        assert period_ordinal(2013, 10, 5, 0,
                              0, 0, 0, 0, get_freq('B')) == 11417
        # Sunday
        assert period_ordinal(2013, 10, 6, 0,
                              0, 0, 0, 0, get_freq('B')) == 11417
        # Monday
        assert period_ordinal(2013, 10, 7, 0,
                              0, 0, 0, 0, get_freq('B')) == 11417
        # Tuesday
        assert period_ordinal(2013, 10, 8, 0,
                              0, 0, 0, 0, get_freq('B')) == 11418
