# Copyright 2010-2023 Kurt McKee <contactme@kurtmckee.org>
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

from .w3dtf import _parse_date_w3dtf

# 8-bit date handling routines written by ytrewq1.
_korean_year = '\ub144' # b3e2 in euc-kr
_korean_month = '\uc6d4' # bff9 in euc-kr
_korean_day = '\uc77c' # c0cf in euc-kr
_korean_am = '\uc624\uc804' # bfc0 c0fc in euc-kr
_korean_pm = '\uc624\ud6c4' # bfc0 c8c4 in euc-kr

_korean_onblog_date_re = re.compile(
    r'(\d{4})%s\s+(\d{2})%s\s+(\d{2})%s\s+(\d{2}):(\d{2}):(\d{2})'
    % (_korean_year, _korean_month, _korean_day)
)

_korean_nate_date_re = re.compile(
    r'(\d{4})-(\d{2})-(\d{2})\s+(%s|%s)\s+(\d{,2}):(\d{,2}):(\d{,2})'
    % (_korean_am, _korean_pm))


def _parse_date_onblog(dateString):
    """Parse a string according to the OnBlog 8-bit date format"""
    m = _korean_onblog_date_re.match(dateString)
    if not m:
        return
    w3dtfdate = '%(year)s-%(month)s-%(day)sT%(hour)s:%(minute)s:%(second)s%(zonediff)s' % \
                {'year': m.group(1), 'month': m.group(2), 'day': m.group(3),
                 'hour': m.group(4), 'minute': m.group(5), 'second': m.group(6),
                 'zonediff': '+09:00'}
    return _parse_date_w3dtf(w3dtfdate)


def _parse_date_nate(dateString):
    """Parse a string according to the Nate 8-bit date format"""
    m = _korean_nate_date_re.match(dateString)
    if not m:
        return
    hour = int(m.group(5))
    ampm = m.group(4)
    if ampm == _korean_pm:
        hour += 12
    hour = str(hour)
    if len(hour) == 1:
        hour = '0' + hour
    w3dtfdate = '%(year)s-%(month)s-%(day)sT%(hour)s:%(minute)s:%(second)s%(zonediff)s' % \
                {
                    'year': m.group(1),
                    'month': m.group(2),
                    'day': m.group(3),
                    'hour': hour,
                    'minute': m.group(6),
                    'second': m.group(7),
                    'zonediff': '+09:00',
                }
    return _parse_date_w3dtf(w3dtfdate)
