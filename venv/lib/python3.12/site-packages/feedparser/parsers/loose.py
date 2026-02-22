# The loose feed parser that interfaces with an SGML parsing library
# Copyright 2010-2025 Kurt McKee <contactme@kurtmckee.org>
# Copyright 2002-2008 Mark Pilgrim
# All rights reserved.
#
# This file is a part of feedparser.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
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

class _LooseFeedParser(object):
    contentparams = None

    def __init__(self, baseuri=None, baselang=None, encoding=None, entities=None):
        self.baseuri = baseuri or ''
        self.lang = baselang or None
        self.encoding = encoding or 'utf-8'  # character encoding
        self.entities = entities or {}
        super(_LooseFeedParser, self).__init__()

    @staticmethod
    def _normalize_attributes(kv):
        k = kv[0].lower()
        v = k in ('rel', 'type') and kv[1].lower() or kv[1]
        # the sgml parser doesn't handle entities in attributes, nor
        # does it pass the attribute values through as unicode, while
        # strict xml parsers do -- account for this difference
        v = v.replace('&amp;', '&')
        return k, v

    def decode_entities(self, element, data):
        data = data.replace('&#60;', '&lt;')
        data = data.replace('&#x3c;', '&lt;')
        data = data.replace('&#x3C;', '&lt;')
        data = data.replace('&#62;', '&gt;')
        data = data.replace('&#x3e;', '&gt;')
        data = data.replace('&#x3E;', '&gt;')
        data = data.replace('&#38;', '&amp;')
        data = data.replace('&#x26;', '&amp;')
        data = data.replace('&#34;', '&quot;')
        data = data.replace('&#x22;', '&quot;')
        data = data.replace('&#39;', '&apos;')
        data = data.replace('&#x27;', '&apos;')
        if not self.contentparams.get('type', 'xml').endswith('xml'):
            data = data.replace('&lt;', '<')
            data = data.replace('&gt;', '>')
            data = data.replace('&amp;', '&')
            data = data.replace('&quot;', '"')
            data = data.replace('&apos;', "'")
            data = data.replace('&#x2f;', '/')
            data = data.replace('&#x2F;', '/')
        return data

    @staticmethod
    def strattrs(attrs):
        return ''.join(
            ' %s="%s"' % (n, v.replace('"', '&quot;'))
            for n, v in attrs
        )
