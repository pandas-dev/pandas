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

from .asctime import _parse_date_asctime
from .greek import _parse_date_greek
from .hungarian import _parse_date_hungarian
from .iso8601 import _parse_date_iso8601
from .korean import _parse_date_onblog, _parse_date_nate
from .perforce import _parse_date_perforce
from .rfc822 import _parse_date_rfc822
from .w3dtf import _parse_date_w3dtf

_date_handlers = []


def registerDateHandler(func):
    """Register a date handler function (takes string, returns 9-tuple date in GMT)"""
    _date_handlers.insert(0, func)


def _parse_date(date_string):
    """Parses a variety of date formats into a 9-tuple in GMT"""
    if not date_string:
        return None
    for handler in _date_handlers:
        try:
            date9tuple = handler(date_string)
        except (KeyError, OverflowError, ValueError, AttributeError):
            continue
        if not date9tuple:
            continue
        if len(date9tuple) != 9:
            continue
        return date9tuple
    return None


registerDateHandler(_parse_date_onblog)
registerDateHandler(_parse_date_nate)
registerDateHandler(_parse_date_greek)
registerDateHandler(_parse_date_hungarian)
registerDateHandler(_parse_date_perforce)
registerDateHandler(_parse_date_asctime)
registerDateHandler(_parse_date_iso8601)
registerDateHandler(_parse_date_rfc822)
registerDateHandler(_parse_date_w3dtf)
