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
import urllib.parse

from .html import _BaseHTMLProcessor

# If you want feedparser to allow all URL schemes, set this to ()
# List culled from Python's urlparse documentation at:
#   http://docs.python.org/library/urlparse.html
# as well as from "URI scheme" at Wikipedia:
#   https://secure.wikimedia.org/wikipedia/en/wiki/URI_scheme
# Many more will likely need to be added!
ACCEPTABLE_URI_SCHEMES = (
    'file', 'ftp', 'gopher', 'h323', 'hdl', 'http', 'https', 'imap', 'magnet',
    'mailto', 'mms', 'news', 'nntp', 'prospero', 'rsync', 'rtsp', 'rtspu',
    'sftp', 'shttp', 'sip', 'sips', 'snews', 'svn', 'svn+ssh', 'telnet',
    'wais',
    # Additional common-but-unofficial schemes
    'aim', 'callto', 'cvs', 'facetime', 'feed', 'git', 'gtalk', 'irc', 'ircs',
    'irc6', 'itms', 'mms', 'msnim', 'skype', 'ssh', 'smb', 'svn', 'ymsg',
)

_urifixer = re.compile('^([A-Za-z][A-Za-z0-9+-.]*://)(/*)(.*?)')


def _urljoin(base, uri):
    uri = _urifixer.sub(r'\1\3', uri)
    try:
        uri = urllib.parse.urljoin(base, uri)
    except ValueError:
        uri = ''
    return uri


def convert_to_idn(url):
    """Convert a URL to IDN notation"""
    # this function should only be called with a unicode string
    # strategy: if the host cannot be encoded in ascii, then
    # it'll be necessary to encode it in idn form
    parts = list(urllib.parse.urlsplit(url))
    try:
        parts[1].encode('ascii')
    except UnicodeEncodeError:
        # the url needs to be converted to idn notation
        host = parts[1].rsplit(':', 1)
        newhost = []
        port = ''
        if len(host) == 2:
            port = host.pop()
        for h in host[0].split('.'):
            newhost.append(h.encode('idna').decode('utf-8'))
        parts[1] = '.'.join(newhost)
        if port:
            parts[1] += ':' + port
        return urllib.parse.urlunsplit(parts)
    else:
        return url


def make_safe_absolute_uri(base, rel=None):
    # bail if ACCEPTABLE_URI_SCHEMES is empty
    if not ACCEPTABLE_URI_SCHEMES:
        return _urljoin(base, rel or '')
    if not base:
        return rel or ''
    if not rel:
        try:
            scheme = urllib.parse.urlparse(base)[0]
        except ValueError:
            return ''
        if not scheme or scheme in ACCEPTABLE_URI_SCHEMES:
            return base
        return ''
    uri = _urljoin(base, rel)
    if uri.strip().split(':', 1)[0] not in ACCEPTABLE_URI_SCHEMES:
        return ''
    return uri


class RelativeURIResolver(_BaseHTMLProcessor):
    relative_uris = {
        ('a', 'href'),
        ('applet', 'codebase'),
        ('area', 'href'),
        ('audio', 'src'),
        ('blockquote', 'cite'),
        ('body', 'background'),
        ('del', 'cite'),
        ('form', 'action'),
        ('frame', 'longdesc'),
        ('frame', 'src'),
        ('iframe', 'longdesc'),
        ('iframe', 'src'),
        ('head', 'profile'),
        ('img', 'longdesc'),
        ('img', 'src'),
        ('img', 'usemap'),
        ('input', 'src'),
        ('input', 'usemap'),
        ('ins', 'cite'),
        ('link', 'href'),
        ('object', 'classid'),
        ('object', 'codebase'),
        ('object', 'data'),
        ('object', 'usemap'),
        ('q', 'cite'),
        ('script', 'src'),
        ('source', 'src'),
        ('video', 'poster'),
        ('video', 'src'),
    }

    def __init__(self, baseuri, encoding, _type):
        _BaseHTMLProcessor.__init__(self, encoding, _type)
        self.baseuri = baseuri

    def resolve_uri(self, uri):
        return make_safe_absolute_uri(self.baseuri, uri.strip())

    def unknown_starttag(self, tag, attrs):
        attrs = self.normalize_attrs(attrs)
        attrs = [(key, ((tag, key) in self.relative_uris) and self.resolve_uri(value) or value) for key, value in attrs]
        super(RelativeURIResolver, self).unknown_starttag(tag, attrs)


def resolve_relative_uris(html_source, base_uri, encoding, type_):
    p = RelativeURIResolver(base_uri, encoding, type_)
    p.feed(html_source)
    return p.output()
