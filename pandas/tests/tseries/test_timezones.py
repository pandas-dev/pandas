# pylint: disable-msg=E1101,W0612
import pytest

import pytz
import dateutil

from datetime import datetime

from pandas._libs.tslibs import timezones
from pandas import Timestamp


class TestTimeZoneSupportPytz(object):

    def tz(self, tz):
        # Construct a timezone object from a string. Overridden in subclass to
        # parameterize tests.
        return pytz.timezone(tz)

    def tzstr(self, tz):
        # Construct a timezone string from a string. Overridden in subclass to
        # parameterize tests.
        return tz

    def localize(self, tz, x):
        return tz.localize(x)

    def normalize(self, ts):
        tzinfo = ts.tzinfo
        return tzinfo.normalize(ts)

    def cmptz(self, tz1, tz2):
        # Compare two timezones. Overridden in subclass to parameterize
        # tests.
        return tz1.zone == tz2.zone

    # test utility methods
    def test_infer_tz(self):
        eastern = self.tz('US/Eastern')
        utc = pytz.utc

        _start = datetime(2001, 1, 1)
        _end = datetime(2009, 1, 1)

        start = self.localize(eastern, _start)
        end = self.localize(eastern, _end)
        assert (timezones.infer_tzinfo(start, end) is
                self.localize(eastern, _start).tzinfo)
        assert (timezones.infer_tzinfo(start, None) is
                self.localize(eastern, _start).tzinfo)
        assert (timezones.infer_tzinfo(None, end) is
                self.localize(eastern, _end).tzinfo)

        start = utc.localize(_start)
        end = utc.localize(_end)
        assert (timezones.infer_tzinfo(start, end) is utc)

        end = self.localize(eastern, _end)
        pytest.raises(Exception, timezones.infer_tzinfo, start, end)
        pytest.raises(Exception, timezones.infer_tzinfo, end, start)

    def test_replace_across_dst(self):
        # GH#18319 check that 1) timezone is correctly normalized and
        # 2) that hour is not incorrectly changed by this normalization
        tz = self.tz('US/Eastern')

        ts_naive = Timestamp('2017-12-03 16:03:30')
        ts_aware = self.localize(tz, ts_naive)

        # Preliminary sanity-check
        assert ts_aware == self.normalize(ts_aware)

        # Replace across DST boundary
        ts2 = ts_aware.replace(month=6)

        # Check that `replace` preserves hour literal
        assert (ts2.hour, ts2.minute) == (ts_aware.hour, ts_aware.minute)

        # Check that post-replace object is appropriately normalized
        ts2b = self.normalize(ts2)
        assert ts2 == ts2b


class TestTimeZoneSupportDateutil(TestTimeZoneSupportPytz):

    def tz(self, tz):
        """
        Construct a dateutil timezone.
        Use tslib.maybe_get_tz so that we get the filename on the tz right
        on windows. See #7337.
        """
        return timezones.maybe_get_tz('dateutil/' + tz)

    def tzstr(self, tz):
        """ Construct a timezone string from a string. Overridden in subclass
        to parameterize tests. """
        return 'dateutil/' + tz

    def cmptz(self, tz1, tz2):
        """ Compare two timezones. Overridden in subclass to parameterize
        tests. """
        return tz1 == tz2

    def localize(self, tz, x):
        return x.replace(tzinfo=tz)

    def normalize(self, ts):
        # no-op for dateutil
        return ts
