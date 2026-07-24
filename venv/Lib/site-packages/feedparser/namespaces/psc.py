# Support for the Podlove Simple Chapters format
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

import datetime
import re

from .. import util


class Namespace(object):
    supported_namespaces = {
        'http://podlove.org/simple-chapters': 'psc',
    }

    def __init__(self):
        # chapters will only be captured while psc_chapters_flag is True.
        self.psc_chapters_flag = False
        super(Namespace, self).__init__()

    def _start_psc_chapters(self, attrs_d):
        context = self._get_context()
        if 'psc_chapters' not in context:
            self.psc_chapters_flag = True
            attrs_d['chapters'] = []
            context['psc_chapters'] = util.FeedParserDict(attrs_d)

    def _end_psc_chapters(self):
        self.psc_chapters_flag = False

    def _start_psc_chapter(self, attrs_d):
        if self.psc_chapters_flag:
            start = self._get_attribute(attrs_d, 'start')
            attrs_d['start_parsed'] = _parse_psc_chapter_start(start)

            context = self._get_context()['psc_chapters']
            context['chapters'].append(util.FeedParserDict(attrs_d))


format_ = re.compile(r'^((\d{2}):)?(\d{2}):(\d{2})(\.(\d{3}))?$')


def _parse_psc_chapter_start(start):
    m = format_.match(start)
    if m is None:
        return None

    _, h, m, s, _, ms = m.groups()
    h, m, s, ms = (int(h or 0), int(m), int(s), int(ms or 0))
    return datetime.timedelta(0, h*60*60 + m*60 + s, ms*1000)
