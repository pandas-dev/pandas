# Copyright 2010-2025 Kurt McKee <contactme@kurtmckee.org>
# Copyright 2002-2008 Mark Pilgrim
# All rights reserved.
#
# This file is a part of feedparser.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import re
import time

# ISO-8601 date parsing routines written by Fazal Majid.
# The ISO 8601 standard is very convoluted and irregular - a full ISO 8601
# parser is beyond the scope of feedparser and would be a worthwhile addition
# to the Python library.
# A single regular expression cannot parse ISO 8601 date formats into groups
# as the standard is highly irregular (for instance is 030104 2003-01-04 or
# 0301-04-01), so we use templates instead.
# Please note the order in templates is significant because we need a
# greedy match.
_iso8601_tmpl = [
    'YYYY-?MM-?DD',
    'YYYY-0MM?-?DD',
    'YYYY-MM',
    'YYYY-?OOO',
    'YY-?MM-?DD',
    'YY-?OOO',
    'YYYY',
    '-YY-?MM',
    '-OOO',
    '-YY',
    '--MM-?DD',
    '--MM',
    '---DD',
    'CC',
    '',
]

_iso8601_re = [
    tmpl.replace(
    'YYYY', r'(?P<year>\d{4})').replace(
    'YY', r'(?P<year>\d\d)').replace(
    'MM', r'(?P<month>[01]\d)').replace(
    'DD', r'(?P<day>[0123]\d)').replace(
    'OOO', r'(?P<ordinal>[0123]\d\d)').replace(
    'CC', r'(?P<century>\d\d$)')
    + r'(T?(?P<hour>\d{2}):(?P<minute>\d{2})'
    + r'(:(?P<second>\d{2}))?'
    + r'(\.(?P<fracsecond>\d+))?'
    + r'(?P<tz>[+-](?P<tzhour>\d{2})(:(?P<tzmin>\d{2}))?|Z)?)?'
    for tmpl in _iso8601_tmpl]
try:
    del tmpl
except NameError:
    pass
_iso8601_matches = [re.compile(regex).match for regex in _iso8601_re]
try:
    del regex
except NameError:
    pass


def _parse_date_iso8601(date_string):
    """Parse a variety of ISO-8601-compatible formats like 20040105"""
    m = None
    for _iso8601_match in _iso8601_matches:
        m = _iso8601_match(date_string)
        if m:
            break
    if not m:
        return
    if m.span() == (0, 0):
        return
    params = m.groupdict()
    ordinal = params.get('ordinal', 0)
    if ordinal:
        ordinal = int(ordinal)
    else:
        ordinal = 0
    year = params.get('year', '--')
    if not year or year == '--':
        year = time.gmtime()[0]
    elif len(year) == 2:
        # ISO 8601 assumes current century, i.e. 93 -> 2093, NOT 1993
        year = 100 * int(time.gmtime()[0] / 100) + int(year)
    else:
        year = int(year)
    month = params.get('month', '-')
    if not month or month == '-':
        # ordinals are NOT normalized by mktime, we simulate them
        # by setting month=1, day=ordinal
        if ordinal:
            month = 1
        else:
            month = time.gmtime()[1]
    month = int(month)
    day = params.get('day', 0)
    if not day:
        # see above
        if ordinal:
            day = ordinal
        elif params.get('century', 0) or \
                 params.get('year', 0) or params.get('month', 0):
            day = 1
        else:
            day = time.gmtime()[2]
    else:
        day = int(day)
    # special case of the century - is the first year of the 21st century
    # 2000 or 2001 ? The debate goes on...
    if 'century' in params:
        year = (int(params['century']) - 1) * 100 + 1
    # in ISO 8601 most fields are optional
    for field in ['hour', 'minute', 'second', 'tzhour', 'tzmin']:
        if not params.get(field, None):
            params[field] = 0
    hour = int(params.get('hour', 0))
    minute = int(params.get('minute', 0))
    second = int(float(params.get('second', 0)))
    # weekday is normalized by mktime(), we can ignore it
    weekday = 0
    daylight_savings_flag = -1
    tm = [year, month, day, hour, minute, second, weekday,
          ordinal, daylight_savings_flag]
    # ISO 8601 time zone adjustments
    tz = params.get('tz')
    if tz and tz != 'Z':
        if tz[0] == '-':
            tm[3] += int(params.get('tzhour', 0))
            tm[4] += int(params.get('tzmin', 0))
        elif tz[0] == '+':
            tm[3] -= int(params.get('tzhour', 0))
            tm[4] -= int(params.get('tzmin', 0))
        else:
            return None
    # Python's time.mktime() is a wrapper around the ANSI C mktime(3c)
    # which is guaranteed to normalize d/m/y/h/m/s.
    # Many implementations have bugs, but we'll pretend they don't.
    return time.localtime(time.mktime(tuple(tm)))
