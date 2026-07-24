# Support for the Creative Commons licensing extensions
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

from ..util import FeedParserDict


class Namespace(object):
    supported_namespaces = {
        # RDF-based namespace
        'http://creativecommons.org/ns#license': 'cc',

        # Old RDF-based namespace
        'http://web.resource.org/cc/': 'cc',

        # RSS-based namespace
        'http://cyber.law.harvard.edu/rss/creativeCommonsRssModule.html': 'creativecommons',

        # Old RSS-based namespace
        'http://backend.userland.com/creativeCommonsRssModule': 'creativecommons',
    }

    def _start_cc_license(self, attrs_d):
        context = self._get_context()
        value = self._get_attribute(attrs_d, 'rdf:resource')
        attrs_d = FeedParserDict()
        attrs_d['rel'] = 'license'
        if value:
            attrs_d['href'] = value
        context.setdefault('links', []).append(attrs_d)

    def _start_creativecommons_license(self, attrs_d):
        self.push('license', 1)
    _start_creativeCommons_license = _start_creativecommons_license

    def _end_creativecommons_license(self):
        value = self.pop('license')
        context = self._get_context()
        attrs_d = FeedParserDict()
        attrs_d['rel'] = 'license'
        if value:
            attrs_d['href'] = value
        context.setdefault('links', []).append(attrs_d)
        del context['license']
    _end_creativeCommons_license = _end_creativecommons_license
