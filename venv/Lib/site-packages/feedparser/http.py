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

import base64
import datetime
import gzip
import io
import re
import struct
import urllib.parse
import urllib.request
import zlib

from .datetimes import _parse_date
from .urls import convert_to_idn


# HTTP "Accept" header to send to servers when downloading feeds.  If you don't
# want to send an Accept header, set this to None.
ACCEPT_HEADER = "application/atom+xml,application/rdf+xml,application/rss+xml,application/x-netcdf,application/xml;q=0.9,text/xml;q=0.2,*/*;q=0.1"


class _FeedURLHandler(urllib.request.HTTPDigestAuthHandler, urllib.request.HTTPRedirectHandler, urllib.request.HTTPDefaultErrorHandler):
    def http_error_default(self, req, fp, code, msg, headers):
        # The default implementation just raises HTTPError.
        # Forget that.
        fp.status = code
        return fp

    def http_error_301(self, req, fp, code, msg, hdrs):
        result = urllib.request.HTTPRedirectHandler.http_error_301(self, req, fp, code, msg, hdrs)
        if not result:
            return fp
        result.status = code
        result.newurl = result.geturl()
        return result

    # The default implementations in urllib.request.HTTPRedirectHandler
    # are identical, so hardcoding a http_error_301 call above
    # won't affect anything
    http_error_300 = http_error_301
    http_error_302 = http_error_301
    http_error_303 = http_error_301
    http_error_307 = http_error_301

    def http_error_401(self, req, fp, code, msg, headers):
        # Check if
        # - server requires digest auth, AND
        # - we tried (unsuccessfully) with basic auth, AND
        # If all conditions hold, parse authentication information
        # out of the Authorization header we sent the first time
        # (for the username and password) and the WWW-Authenticate
        # header the server sent back (for the realm) and retry
        # the request with the appropriate digest auth headers instead.
        # This evil genius hack has been brought to you by Aaron Swartz.
        host = urllib.parse.urlparse(req.get_full_url())[1]
        if 'Authorization' not in req.headers or 'WWW-Authenticate' not in headers:
            return self.http_error_default(req, fp, code, msg, headers)
        auth = base64.decodebytes(req.headers['Authorization'].split(' ')[1].encode()).decode()
        user, passw = auth.split(':')
        realm = re.findall('realm="([^"]*)"', headers['WWW-Authenticate'])[0]
        self.add_password(realm, host, user, passw)
        retry = self.http_error_auth_reqed('www-authenticate', host, req, headers)
        self.reset_retry_count()
        return retry


def _build_urllib2_request(url, agent, accept_header, etag, modified, referrer, auth, request_headers):
    request = urllib.request.Request(url)
    request.add_header('User-Agent', agent)
    if etag:
        request.add_header('If-None-Match', etag)
    if isinstance(modified, str):
        modified = _parse_date(modified)
    elif isinstance(modified, datetime.datetime):
        modified = modified.utctimetuple()
    if modified:
        # format into an RFC 1123-compliant timestamp. We can't use
        # time.strftime() since the %a and %b directives can be affected
        # by the current locale, but RFC 2616 states that dates must be
        # in English.
        short_weekdays = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']
        months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        request.add_header('If-Modified-Since', '%s, %02d %s %04d %02d:%02d:%02d GMT' % (short_weekdays[modified[6]], modified[2], months[modified[1] - 1], modified[0], modified[3], modified[4], modified[5]))
    if referrer:
        request.add_header('Referer', referrer)
    request.add_header('Accept-encoding', 'gzip, deflate')
    if auth:
        request.add_header('Authorization', 'Basic %s' % auth)
    if accept_header:
        request.add_header('Accept', accept_header)
    # use this for whatever -- cookies, special headers, etc
    # [('Cookie','Something'),('x-special-header','Another Value')]
    for header_name, header_value in request_headers.items():
        request.add_header(header_name, header_value)
    request.add_header('A-IM', 'feed')  # RFC 3229 support
    return request


def get(url, etag=None, modified=None, agent=None, referrer=None, handlers=None, request_headers=None, result=None):
    if handlers is None:
        handlers = []
    elif not isinstance(handlers, list):
        handlers = [handlers]
    if request_headers is None:
        request_headers = {}

    # Deal with the feed URI scheme
    if url.startswith('feed:http'):
        url = url[5:]
    elif url.startswith('feed:'):
        url = 'http:' + url[5:]
    if not agent:
        from . import USER_AGENT
        agent = USER_AGENT
    # Test for inline user:password credentials for HTTP basic auth
    auth = None
    if not url.startswith('ftp:'):
        url_pieces = urllib.parse.urlparse(url)
        if url_pieces.username:
            new_pieces = list(url_pieces)
            new_pieces[1] = url_pieces.hostname
            if url_pieces.port:
                new_pieces[1] = f'{url_pieces.hostname}:{url_pieces.port}'
            url = urllib.parse.urlunparse(new_pieces)
            auth = base64.standard_b64encode(f'{url_pieces.username}:{url_pieces.password}'.encode()).decode()

    # iri support
    if not isinstance(url, bytes):
        url = convert_to_idn(url)

    # Prevent UnicodeEncodeErrors caused by Unicode characters in the path.
    bits = []
    for c in url:
        try:
            c.encode('ascii')
        except UnicodeEncodeError:
            bits.append(urllib.parse.quote(c))
        else:
            bits.append(c)
    url = ''.join(bits)

    # try to open with urllib2 (to use optional headers)
    request = _build_urllib2_request(url, agent, ACCEPT_HEADER, etag, modified, referrer, auth, request_headers)
    opener = urllib.request.build_opener(*tuple(handlers + [_FeedURLHandler()]))
    opener.addheaders = []  # RMK - must clear so we only send our custom User-Agent
    f = opener.open(request)
    data = f.read()
    f.close()

    # lowercase all of the HTTP headers for comparisons per RFC 2616
    result['headers'] = {k.lower(): v for k, v in f.headers.items()}

    # if feed is gzip-compressed, decompress it
    if data and 'gzip' in result['headers'].get('content-encoding', ''):
        try:
            data = gzip.GzipFile(fileobj=io.BytesIO(data)).read()
        except (EOFError, IOError, struct.error) as e:
            # IOError can occur if the gzip header is bad.
            # struct.error can occur if the data is damaged.
            result['bozo'] = True
            result['bozo_exception'] = e
            if isinstance(e, struct.error):
                # A gzip header was found but the data is corrupt.
                # Ideally, we should re-request the feed without the
                # 'Accept-encoding: gzip' header, but we don't.
                data = None
    elif data and 'deflate' in result['headers'].get('content-encoding', ''):
        try:
            data = zlib.decompress(data)
        except zlib.error:
            try:
                # The data may have no headers and no checksum.
                data = zlib.decompress(data, -15)
            except zlib.error as e:
                result['bozo'] = True
                result['bozo_exception'] = e

    # save HTTP headers
    if 'etag' in result['headers']:
        etag = result['headers'].get('etag', '')
        if isinstance(etag, bytes):
            etag = etag.decode('utf-8', 'ignore')
        if etag:
            result['etag'] = etag
    if 'last-modified' in result['headers']:
        modified = result['headers'].get('last-modified', '')
        if modified:
            result['modified'] = modified
            result['modified_parsed'] = _parse_date(modified)
    if isinstance(f.url, bytes):
        result['href'] = f.url.decode('utf-8', 'ignore')
    else:
        result['href'] = f.url
    result['status'] = getattr(f, 'status', None) or 200

    # Stop processing if the server sent HTTP 304 Not Modified.
    if getattr(f, 'code', 0) == 304:
        result['version'] = ''
        result['debug_message'] = 'The feed has not changed since you last checked, ' + \
            'so the server sent no data.  This is a feature, not a bug!'

    return data
