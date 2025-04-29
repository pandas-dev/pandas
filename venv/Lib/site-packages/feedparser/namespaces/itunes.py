# Support for the iTunes format
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
        'http://www.itunes.com/DTDs/PodCast-1.0.dtd': 'itunes',

        # Extra namespace
        'http://example.com/DTDs/PodCast-1.0.dtd': 'itunes',
    }

    def _start_itunes_author(self, attrs_d):
        self._start_author(attrs_d)

    def _end_itunes_author(self):
        self._end_author()

    def _end_itunes_category(self):
        self._end_category()

    def _start_itunes_name(self, attrs_d):
        self._start_name(attrs_d)

    def _end_itunes_name(self):
        self._end_name()

    def _start_itunes_email(self, attrs_d):
        self._start_email(attrs_d)

    def _end_itunes_email(self):
        self._end_email()

    def _start_itunes_subtitle(self, attrs_d):
        self._start_subtitle(attrs_d)

    def _end_itunes_subtitle(self):
        self._end_subtitle()

    def _start_itunes_summary(self, attrs_d):
        self._start_summary(attrs_d)

    def _end_itunes_summary(self):
        self._end_summary()

    def _start_itunes_owner(self, attrs_d):
        self.inpublisher = 1
        self.push('publisher', 0)

    def _end_itunes_owner(self):
        self.pop('publisher')
        self.inpublisher = 0
        self._sync_author_detail('publisher')

    def _end_itunes_keywords(self):
        for term in self.pop('itunes_keywords').split(','):
            if term.strip():
                self._add_tag(term.strip(), 'http://www.itunes.com/', None)

    def _start_itunes_category(self, attrs_d):
        self._add_tag(attrs_d.get('text'), 'http://www.itunes.com/', None)
        self.push('category', 1)

    def _start_itunes_image(self, attrs_d):
        self.push('itunes_image', 0)
        if attrs_d.get('href'):
            self._get_context()['image'] = FeedParserDict({'href': attrs_d.get('href')})
        elif attrs_d.get('url'):
            self._get_context()['image'] = FeedParserDict({'href': attrs_d.get('url')})
    _start_itunes_link = _start_itunes_image

    def _end_itunes_block(self):
        value = self.pop('itunes_block', 0)
        self._get_context()['itunes_block'] = (value == 'yes' or value == 'Yes') and 1 or 0

    def _end_itunes_explicit(self):
        value = self.pop('itunes_explicit', 0)
        # Convert 'yes' -> True, 'clean' to False, and any other value to None
        # False and None both evaluate as False, so the difference can be ignored
        # by applications that only need to know if the content is explicit.
        self._get_context()['itunes_explicit'] = (None, False, True)[(value == 'yes' and 2) or value == 'clean' or 0]
