# -*- coding: utf-8 -*-
import datetime


def get_system_offset():
    """Get system's timezone offset using built-in library time.

    For the Timezone constants (altzone, daylight, timezone, and tzname), the
    value is determined by the timezone rules in effect at module load time or
    the last time tzset() is called and may be incorrect for times in the past.

    To keep compatibility with Windows, we're always importing time module here.
    """
    import time
    if time.daylight and time.localtime().tm_isdst > 0:
        return -time.altzone
    else:
        return -time.timezone


def get_tz_offset(tz):
    """Get timezone's offset using built-in function datetime.utcoffset()."""
    return int(datetime.datetime.now(tz).utcoffset().total_seconds())


def assert_tz_offset(tz):
    """Assert that system's timezone offset equals to the timezone offset found.

    If they don't match, we probably have a misconfiguration, for example, an
    incorrect timezone set in /etc/timezone file in systemd distributions."""
    tz_offset = get_tz_offset(tz)
    system_offset = get_system_offset()
    if tz_offset != system_offset:
        msg = ('Timezone offset does not match system offset: {0} != {1}. '
               'Please, check your config files.').format(
                   tz_offset, system_offset
               )
        raise ValueError(msg)
