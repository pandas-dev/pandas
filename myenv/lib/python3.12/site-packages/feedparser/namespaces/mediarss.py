# Support for the Media RSS format
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

from ..util import FeedParserDict


class Namespace(object):
    supported_namespaces = {
        # Canonical namespace
        'http://search.yahoo.com/mrss/': 'media',

        # Old namespace (no trailing slash)
        'http://search.yahoo.com/mrss': 'media',
    }

    def _start_media_category(self, attrs_d):
        attrs_d.setdefault('scheme', 'http://search.yahoo.com/mrss/category_schema')
        self._start_category(attrs_d)

    def _end_media_category(self):
        self._end_category()

    def _end_media_keywords(self):
        for term in self.pop('media_keywords').split(','):
            if term.strip():
                self._add_tag(term.strip(), None, None)

    def _start_media_title(self, attrs_d):
        self._start_title(attrs_d)

    def _end_media_title(self):
        title_depth = self.title_depth
        self._end_title()
        self.title_depth = title_depth

    def _start_media_group(self, attrs_d):
        # don't do anything, but don't break the enclosed tags either
        pass

    def _start_media_rating(self, attrs_d):
        context = self._get_context()
        context.setdefault('media_rating', attrs_d)
        self.push('rating', 1)

    def _end_media_rating(self):
        rating = self.pop('rating')
        if rating is not None and rating.strip():
            context = self._get_context()
            context['media_rating']['content'] = rating

    def _start_media_credit(self, attrs_d):
        context = self._get_context()
        context.setdefault('media_credit', [])
        context['media_credit'].append(attrs_d)
        self.push('credit', 1)

    def _end_media_credit(self):
        credit = self.pop('credit')
        if credit is not None and credit.strip():
            context = self._get_context()
            context['media_credit'][-1]['content'] = credit

    def _start_media_description(self, attrs_d):
        self._start_description(attrs_d)

    def _end_media_description(self):
        self._end_description()

    def _start_media_restriction(self, attrs_d):
        context = self._get_context()
        context.setdefault('media_restriction', attrs_d)
        self.push('restriction', 1)

    def _end_media_restriction(self):
        restriction = self.pop('restriction')
        if restriction is not None and restriction.strip():
            context = self._get_context()
            context['media_restriction']['content'] = [cc.strip().lower() for cc in restriction.split(' ')]

    def _start_media_license(self, attrs_d):
        context = self._get_context()
        context.setdefault('media_license', attrs_d)
        self.push('license', 1)

    def _end_media_license(self):
        license_ = self.pop('license')
        if license_ is not None and license_.strip():
            context = self._get_context()
            context['media_license']['content'] = license_

    def _start_media_content(self, attrs_d):
        context = self._get_context()
        context.setdefault('media_content', [])
        context['media_content'].append(attrs_d)

    def _start_media_thumbnail(self, attrs_d):
        context = self._get_context()
        context.setdefault('media_thumbnail', [])
        self.push('url', 1) # new
        context['media_thumbnail'].append(attrs_d)

    def _end_media_thumbnail(self):
        url = self.pop('url')
        context = self._get_context()
        if url is not None and url.strip():
            if 'url' not in context['media_thumbnail'][-1]:
                context['media_thumbnail'][-1]['url'] = url

    def _start_media_player(self, attrs_d):
        self.push('media_player', 0)
        self._get_context()['media_player'] = FeedParserDict(attrs_d)

    def _end_media_player(self):
        value = self.pop('media_player')
        context = self._get_context()
        context['media_player']['content'] = value
