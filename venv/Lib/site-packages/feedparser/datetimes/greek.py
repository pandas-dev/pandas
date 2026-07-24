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

from .rfc822 import _parse_date_rfc822

# Unicode strings for Greek date strings
_greek_months = {
   '\u0399\u03b1\u03bd': 'Jan',        # c9e1ed in iso-8859-7
   '\u03a6\u03b5\u03b2': 'Feb',        # d6e5e2 in iso-8859-7
   '\u039c\u03ac\u03ce': 'Mar',        # ccdcfe in iso-8859-7
   '\u039c\u03b1\u03ce': 'Mar',        # cce1fe in iso-8859-7
   '\u0391\u03c0\u03c1': 'Apr',        # c1f0f1 in iso-8859-7
   '\u039c\u03ac\u03b9': 'May',        # ccdce9 in iso-8859-7
   '\u039c\u03b1\u03ca': 'May',        # cce1fa in iso-8859-7
   '\u039c\u03b1\u03b9': 'May',        # cce1e9 in iso-8859-7
   '\u0399\u03bf\u03cd\u03bd': 'Jun',  # c9effded in iso-8859-7
   '\u0399\u03bf\u03bd': 'Jun',        # c9efed in iso-8859-7
   '\u0399\u03bf\u03cd\u03bb': 'Jul',  # c9effdeb in iso-8859-7
   '\u0399\u03bf\u03bb': 'Jul',        # c9f9eb in iso-8859-7
   '\u0391\u03cd\u03b3': 'Aug',        # c1fde3 in iso-8859-7
   '\u0391\u03c5\u03b3': 'Aug',        # c1f5e3 in iso-8859-7
   '\u03a3\u03b5\u03c0': 'Sep',        # d3e5f0 in iso-8859-7
   '\u039f\u03ba\u03c4': 'Oct',        # cfeaf4 in iso-8859-7
   '\u039d\u03bf\u03ad': 'Nov',        # cdefdd in iso-8859-7
   '\u039d\u03bf\u03b5': 'Nov',        # cdefe5 in iso-8859-7
   '\u0394\u03b5\u03ba': 'Dec',        # c4e5ea in iso-8859-7
}

_greek_wdays = {
   '\u039a\u03c5\u03c1': 'Sun',  # caf5f1 in iso-8859-7
   '\u0394\u03b5\u03c5': 'Mon',  # c4e5f5 in iso-8859-7
   '\u03a4\u03c1\u03b9': 'Tue',  # d4f1e9 in iso-8859-7
   '\u03a4\u03b5\u03c4': 'Wed',  # d4e5f4 in iso-8859-7
   '\u03a0\u03b5\u03bc': 'Thu',  # d0e5ec in iso-8859-7
   '\u03a0\u03b1\u03c1': 'Fri',  # d0e1f1 in iso-8859-7
   '\u03a3\u03b1\u03b2': 'Sat',  # d3e1e2 in iso-8859-7
}

_greek_date_format_re = re.compile(r'([^,]+),\s+(\d{2})\s+([^\s]+)\s+(\d{4})\s+(\d{2}):(\d{2}):(\d{2})\s+([^\s]+)')


def _parse_date_greek(date_string):
    """Parse a string according to a Greek 8-bit date format."""
    m = _greek_date_format_re.match(date_string)
    if not m:
        return
    wday = _greek_wdays[m.group(1)]
    month = _greek_months[m.group(3)]
    rfc822date = '%(wday)s, %(day)s %(month)s %(year)s %(hour)s:%(minute)s:%(second)s %(zonediff)s' % \
                 {
                     'wday': wday,
                     'day': m.group(2),
                     'month': month,
                     'year': m.group(4),
                     'hour': m.group(5),
                     'minute': m.group(6),
                     'second': m.group(7),
                     'zonediff': m.group(8),
                 }
    return _parse_date_rfc822(rfc822date)
