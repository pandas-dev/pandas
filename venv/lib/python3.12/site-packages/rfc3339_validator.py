# -*- coding: utf-8 -*-

__author__ = """Nicolas Aimetti"""
__email__ = 'naimetti@yahoo.com.ar'
__version__ = '0.1.4'

import re
import calendar
import six

RFC3339_REGEX_FLAGS = 0
if six.PY3:
    RFC3339_REGEX_FLAGS |= re.ASCII

RFC3339_REGEX = re.compile(r"""
    ^
    (\d{4})      # Year
    -
    (0[1-9]|1[0-2]) # Month
    -
    (\d{2})          # Day
    T
    (?:[01]\d|2[0123]) # Hours
    :
    (?:[0-5]\d)     # Minutes
    :
    (?:[0-5]\d)     # Seconds
    (?:\.\d+)?      # Secfrac
    (?:  Z                              # UTC
       | [+-](?:[01]\d|2[0123]):[0-5]\d # Offset
    )
    $
""", re.VERBOSE | RFC3339_REGEX_FLAGS)


def validate_rfc3339(date_string):
    """
    Validates dates against RFC3339 datetime format
    Leap seconds are no supported.
    """
    m = RFC3339_REGEX.match(date_string)
    if m is None:
        return False
    year, month, day = map(int, m.groups())
    if not year:
        # Year 0 is not valid a valid date
        return False
    (_, max_day) = calendar.monthrange(year, month)
    if not 1 <= day <= max_day:
        return False
    return True
