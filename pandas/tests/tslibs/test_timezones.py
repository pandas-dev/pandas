# -*- coding: utf-8 -*-
from datetime import datetime

import pytest
import pytz
import dateutil.tz

from pandas._libs.tslibs import timezones
from pandas import Timestamp


class TestTimeZoneCacheKey(object):

    @pytest.mark.parametrize('tz_name', list(pytz.common_timezones))
    def test_cache_keys_are_distinct_for_pytz_vs_dateutil(self, tz_name):
        if tz_name == 'UTC':
            # skip utc as it's a special case in dateutil
            return
        tz_p = timezones.maybe_get_tz(tz_name)
        tz_d = timezones.maybe_get_tz('dateutil/' + tz_name)
        if tz_d is None:
            # skip timezones that dateutil doesn't know about.
            return
        assert (timezones._p_tz_cache_key(tz_p) !=
                timezones._p_tz_cache_key(tz_d))


def test_tzlocal():
    # GH#13583
    ts = Timestamp('2011-01-01', tz=dateutil.tz.tzlocal())
    assert ts.tz == dateutil.tz.tzlocal()
    assert "tz='tzlocal()')" in repr(ts)

    tz = timezones.maybe_get_tz('tzlocal()')
    assert tz == dateutil.tz.tzlocal()

    # get offset using normal datetime for test
    offset = dateutil.tz.tzlocal().utcoffset(datetime(2011, 1, 1))
    offset = offset.total_seconds() * 1000000000
    assert ts.value + offset == Timestamp('2011-01-01').value
