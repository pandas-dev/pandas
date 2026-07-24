# Based on https://github.com/MagicStack/httptools
#

from cpython cimport (
    Py_buffer,
    PyBUF_SIMPLE,
    PyBuffer_Release,
    PyBytes_AsString,
    PyBytes_AsStringAndSize,
    PyObject_GetBuffer,
)
from cpython.mem cimport PyMem_Free, PyMem_Malloc
from libc.limits cimport ULLONG_MAX
from libc.string cimport memcpy

from multidict import CIMultiDict as _CIMultiDict, CIMultiDictProxy as _CIMultiDictProxy
from yarl import URL as _URL

from aiohttp import hdrs
from aiohttp.helpers import DEBUG, set_exception

from .http_exceptions import (
    BadHttpMessage,
    BadHttpMethod,
    BadStatusLine,
    ContentLengthError,
    InvalidHeader,
    InvalidURLError,
    LineTooLong,
    PayloadEncodingError,
    TransferEncodingError,
)
from .http_parser import DeflateBuffer as _DeflateBuffer
from .http_writer import (
    HttpVersion as _HttpVersion,
    HttpVersion10 as _HttpVersion10,
    HttpVersion11 as _HttpVersion11,
)
from .streams import EMPTY_PAYLOAD as _EMPTY_PAYLOAD, StreamReader as _StreamReader

cimport cython

from aiohttp cimport _cparser as cparser

include "_headers.pxi"

from aiohttp cimport _find_header


cdef frozenset ALLOWED_UPGRADES = frozenset({"websocket"})
DEF DEFAULT_FREELIST_SIZE = 250

cdef extern from "Python.h":
    int PyByteArray_Resize(object, Py_ssize_t) except -1
    Py_ssize_t PyByteArray_Size(object) except -1
    char* PyByteArray_AsString(object)

__all__ = ('HttpRequestParser', 'HttpResponseParser',
           'RawRequestMessage', 'RawResponseMessage')

cdef object URL = _URL
cdef object URL_build = URL.build
cdef object CIMultiDict = _CIMultiDict
cdef object CIMultiDictProxy = _CIMultiDictProxy
cdef object HttpVersion = _HttpVersion
cdef object HttpVersion10 = _HttpVersion10
cdef object HttpVersion11 = _HttpVersion11
cdef object SEC_WEBSOCKET_KEY1 = hdrs.SEC_WEBSOCKET_KEY1
cdef object CONTENT_ENCODING = hdrs.CONTENT_ENCODING
cdef object EMPTY_PAYLOAD = _EMPTY_PAYLOAD
cdef object StreamReader = _StreamReader
cdef object DeflateBuffer = _DeflateBuffer
cdef tuple EMPTY_FEED_DATA_RESULT = ((), False, b"")

# RFC 9110 singleton headers — duplicates are rejected in strict mode.
# In lax mode (response parser default), the check is skipped entirely
# since real-world servers (e.g. Google APIs, Werkzeug) commonly send
# duplicate headers like Content-Type or Server.
cdef frozenset SINGLETON_HEADERS = frozenset({
    hdrs.CONTENT_LENGTH,
    hdrs.CONTENT_LOCATION,
    hdrs.CONTENT_RANGE,
    hdrs.CONTENT_TYPE,
    hdrs.ETAG,
    hdrs.HOST,
    hdrs.MAX_FORWARDS,
    hdrs.SERVER,
    hdrs.TRANSFER_ENCODING,
    hdrs.USER_AGENT,
})

cdef inline object extend(object buf, const char* at, size_t length):
    cdef Py_ssize_t s
    cdef char* ptr
    s = PyByteArray_Size(buf)
    PyByteArray_Resize(buf, s + length)
    ptr = PyByteArray_AsString(buf)
    memcpy(ptr + s, at, length)


# The method-name table and its length come straight from llhttp's canonical
# HTTP_ALL_METHOD_MAP, so they track the vendored llhttp version automatically
# instead of relying on a hand-maintained method count.
cdef extern from *:
    """
    #include "llhttp.h"

    #define _AIOHTTP_METHOD_NAME(NUM, NAME, STRING) [NUM] = #STRING,
    static const char* const _aiohttp_method_names[] = {
        HTTP_ALL_METHOD_MAP(_AIOHTTP_METHOD_NAME)
    };
    #undef _AIOHTTP_METHOD_NAME
    """
    const char* _aiohttp_method_names[]
    const int METHODS_COUNT "((int)(sizeof(_aiohttp_method_names) / sizeof(_aiohttp_method_names[0])))"


cdef list _http_method = []

for i in range(METHODS_COUNT):
    assert _aiohttp_method_names[i] is not NULL
    _http_method.append(_aiohttp_method_names[i].decode('ascii'))


cdef inline object find_header(bytes raw_header):
    cdef Py_ssize_t size
    cdef char *buf
    cdef int idx
    PyBytes_AsStringAndSize(raw_header, &buf, &size)
    idx = _find_header.find_header(buf, size)
    if idx == -1:
        return raw_header.decode('utf-8', 'surrogateescape')
    return headers[idx]


@cython.freelist(DEFAULT_FREELIST_SIZE)
cdef class RawRequestMessage:
    cdef readonly str method
    cdef readonly str path
    cdef readonly object version  # HttpVersion
    cdef readonly object headers  # CIMultiDict
    cdef readonly object raw_headers  # tuple
    cdef readonly object should_close
    cdef readonly object compression
    cdef readonly object upgrade
    cdef readonly object chunked
    cdef readonly object url  # yarl.URL

    def __init__(self, method, path, version, headers, raw_headers,
                 should_close, compression, upgrade, chunked, url):
        self.method = method
        self.path = path
        self.version = version
        self.headers = headers
        self.raw_headers = raw_headers
        self.should_close = should_close
        self.compression = compression
        self.upgrade = upgrade
        self.chunked = chunked
        self.url = url

    def __repr__(self):
        info = []
        info.append(("method", self.method))
        info.append(("path", self.path))
        info.append(("version", self.version))
        info.append(("headers", self.headers))
        info.append(("raw_headers", self.raw_headers))
        info.append(("should_close", self.should_close))
        info.append(("compression", self.compression))
        info.append(("upgrade", self.upgrade))
        info.append(("chunked", self.chunked))
        info.append(("url", self.url))
        sinfo = ', '.join(name + '=' + repr(val) for name, val in info)
        return '<RawRequestMessage(' + sinfo + ')>'

    def _replace(self, **dct):
        cdef RawRequestMessage ret
        ret = _new_request_message(self.method,
                                   self.path,
                                   self.version,
                                   self.headers,
                                   self.raw_headers,
                                   self.should_close,
                                   self.compression,
                                   self.upgrade,
                                   self.chunked,
                                   self.url)
        if "method" in dct:
            ret.method = dct["method"]
        if "path" in dct:
            ret.path = dct["path"]
        if "version" in dct:
            ret.version = dct["version"]
        if "headers" in dct:
            ret.headers = dct["headers"]
        if "raw_headers" in dct:
            ret.raw_headers = dct["raw_headers"]
        if "should_close" in dct:
            ret.should_close = dct["should_close"]
        if "compression" in dct:
            ret.compression = dct["compression"]
        if "upgrade" in dct:
            ret.upgrade = dct["upgrade"]
        if "chunked" in dct:
            ret.chunked = dct["chunked"]
        if "url" in dct:
            ret.url = dct["url"]
        return ret

cdef _new_request_message(str method,
                           str path,
                           object version,
                           object headers,
                           object raw_headers,
                           bint should_close,
                           object compression,
                           bint upgrade,
                           bint chunked,
                           object url):
    cdef RawRequestMessage ret
    ret = RawRequestMessage.__new__(RawRequestMessage)
    ret.method = method
    ret.path = path
    ret.version = version
    ret.headers = headers
    ret.raw_headers = raw_headers
    ret.should_close = should_close
    ret.compression = compression
    ret.upgrade = upgrade
    ret.chunked = chunked
    ret.url = url
    return ret


@cython.freelist(DEFAULT_FREELIST_SIZE)
cdef class RawResponseMessage:
    cdef readonly object version  # HttpVersion
    cdef readonly int code
    cdef readonly str reason
    cdef readonly object headers  # CIMultiDict
    cdef readonly object raw_headers  # tuple
    cdef readonly object should_close
    cdef readonly object compression
    cdef readonly object upgrade
    cdef readonly object chunked

    def __init__(self, version, code, reason, headers, raw_headers,
                 should_close, compression, upgrade, chunked):
        self.version = version
        self.code = code
        self.reason = reason
        self.headers = headers
        self.raw_headers = raw_headers
        self.should_close = should_close
        self.compression = compression
        self.upgrade = upgrade
        self.chunked = chunked

    def __repr__(self):
        info = []
        info.append(("version", self.version))
        info.append(("code", self.code))
        info.append(("reason", self.reason))
        info.append(("headers", self.headers))
        info.append(("raw_headers", self.raw_headers))
        info.append(("should_close", self.should_close))
        info.append(("compression", self.compression))
        info.append(("upgrade", self.upgrade))
        info.append(("chunked", self.chunked))
        sinfo = ', '.join(name + '=' + repr(val) for name, val in info)
        return '<RawResponseMessage(' + sinfo + ')>'


cdef _new_response_message(object version,
                           int code,
                           str reason,
                           object headers,
                           object raw_headers,
                           bint should_close,
                           object compression,
                           bint upgrade,
                           bint chunked):
    cdef RawResponseMessage ret
    ret = RawResponseMessage.__new__(RawResponseMessage)
    ret.version = version
    ret.code = code
    ret.reason = reason
    ret.headers = headers
    ret.raw_headers = raw_headers
    ret.should_close = should_close
    ret.compression = compression
    ret.upgrade = upgrade
    ret.chunked = chunked
    return ret


@cython.internal
cdef class HttpParser:

    cdef:
        cparser.llhttp_t* _cparser
        cparser.llhttp_settings_t* _csettings

        bytes _raw_name
        object _name
        bytes _raw_value
        bint      _has_value
        int _header_name_size

        readonly object protocol
        object _loop
        object _timer

        size_t _max_line_size
        size_t _max_field_size
        size_t _max_headers
        bint _response_with_body
        bint _read_until_eof
        bint _lax

        bytes   _tail
        bint    _started
        object  _url
        bytearray   _buf
        str     _path
        str     _reason
        list    _headers
        set     _seen_singletons
        list    _raw_headers
        bint    _upgraded
        bint    _pending_upgrade
        list    _messages
        bint    _more_data_available
        bint    _paused
        Py_ssize_t _msg_in_flight
        Py_ssize_t _max_msg_queue_size
        bint    _eof_pending
        object  _payload
        unsigned long long _content_length_expected
        bint    _payload_error
        object  _payload_exception
        object  _last_error
        bint    _auto_decompress
        int     _limit

        str     _content_encoding

        Py_buffer py_buf

    def __cinit__(self):
        self._cparser = <cparser.llhttp_t*> \
                                PyMem_Malloc(sizeof(cparser.llhttp_t))
        if self._cparser is NULL:
            raise MemoryError()

        self._csettings = <cparser.llhttp_settings_t*> \
                                PyMem_Malloc(sizeof(cparser.llhttp_settings_t))
        if self._csettings is NULL:
            raise MemoryError()

    def __dealloc__(self):
        PyMem_Free(self._cparser)
        PyMem_Free(self._csettings)

    cdef _init(
        self, cparser.llhttp_type mode,
        object protocol, object loop, int limit,
        object timer=None,
        size_t max_line_size=8190, size_t max_headers=128,
        size_t max_field_size=8190, payload_exception=None,
        bint response_with_body=True, bint read_until_eof=False,
        bint auto_decompress=True,
        Py_ssize_t max_msg_queue_size=0,
    ):
        cparser.llhttp_settings_init(self._csettings)
        cparser.llhttp_init(self._cparser, mode, self._csettings)
        self._cparser.data = <void*>self
        self._cparser.content_length = 0
        self._content_length_expected = 0

        self.protocol = protocol
        self._loop = loop
        self._timer = timer

        self._buf = bytearray()
        self._more_data_available = False
        self._paused = False
        self._msg_in_flight = 0
        self._max_msg_queue_size = max_msg_queue_size
        self._eof_pending = False
        self._payload = None
        self._payload_error = 0
        self._payload_exception = payload_exception
        self._messages = []

        self._raw_name = b""
        self._raw_value = b""
        self._tail = b""
        self._has_value = False
        self._header_name_size = 0

        self._max_line_size = max_line_size
        self._max_headers = max_headers
        self._max_field_size = max_field_size
        self._response_with_body = response_with_body
        self._read_until_eof = read_until_eof
        self._upgraded = False
        self._pending_upgrade = False
        self._auto_decompress = auto_decompress
        self._content_encoding = None
        self._lax = False
        self._seen_singletons = set()

        self._csettings.on_url = cb_on_url
        self._csettings.on_status = cb_on_status
        self._csettings.on_header_field = cb_on_header_field
        self._csettings.on_header_value = cb_on_header_value
        self._csettings.on_headers_complete = cb_on_headers_complete
        self._csettings.on_body = cb_on_body
        self._csettings.on_message_begin = cb_on_message_begin
        self._csettings.on_message_complete = cb_on_message_complete
        self._csettings.on_chunk_header = cb_on_chunk_header
        self._csettings.on_chunk_complete = cb_on_chunk_complete

        self._last_error = None
        self._limit = limit

    cdef _process_header(self):
        cdef str value
        if self._raw_name != b"":
            name = find_header(self._raw_name)
            value = self._raw_value.decode('utf-8', 'surrogateescape')

            # reject null bytes in header values - matches the Python parser
            # check at http_parser.py. llhttp in lenient mode doesn't reject
            # these itself, so we need to catch them here.
            # ref: RFC 9110 section 5.5 (CTL chars forbidden in field values)
            if "\x00" in value:
                raise InvalidHeader(self._raw_value)

            if not self._lax and name in SINGLETON_HEADERS:
                if name in self._seen_singletons:
                    raise BadHttpMessage(f"Duplicate '{name}' header found.")
                self._seen_singletons.add(name)
            self._headers.append((name, value))
            if len(self._headers) > self._max_headers:
                raise BadHttpMessage("Too many headers received")

            if name is CONTENT_ENCODING:
                self._content_encoding = value

            self._has_value = False
            self._header_name_size = 0
            self._raw_headers.append((self._raw_name, self._raw_value))
            self._raw_name = b""
            self._raw_value = b""

    cdef _on_header_field(self, char* at, size_t length):
        if self._has_value:
            self._process_header()

        if self._raw_name == b"":
            self._raw_name = at[:length]
        else:
            self._raw_name += at[:length]

    cdef _on_header_value(self, char* at, size_t length):
        if self._raw_value == b"":
            self._raw_value = at[:length]
        else:
            self._raw_value += at[:length]
        self._has_value = True

    cdef _on_headers_complete(self):
        cdef str h_upg
        cdef str enc

        self._process_header()

        http_version = self.http_version()
        should_close = not cparser.llhttp_should_keep_alive(self._cparser)
        upgrade = self._cparser.upgrade
        chunked = self._cparser.flags & cparser.F_CHUNKED

        raw_headers = tuple(self._raw_headers)
        headers = CIMultiDictProxy(CIMultiDict(self._headers))

        if self._cparser.type == cparser.HTTP_REQUEST:
            if http_version == HttpVersion11 and hdrs.HOST not in headers:
                raise BadHttpMessage("Missing 'Host' header in request.")
            h_upg = headers.get("upgrade", "")
            if (upgrade and h_upg.isascii() and h_upg.lower() in ALLOWED_UPGRADES) or self._cparser.method == cparser.HTTP_CONNECT:
                # https://www.rfc-editor.org/info/rfc9110/#section-7.8-15
                # Defer the protocol switch until the complete request has been
                # received.
                self._pending_upgrade = True
        else:
            if upgrade and self._cparser.status_code == 101:
                # llhttp pauses for a 101 on its own; just mark the pending
                # switch so feed_data returns the upgraded-protocol tail.
                self._pending_upgrade = True

        # do not support old websocket spec
        if SEC_WEBSOCKET_KEY1 in headers:
            raise InvalidHeader(SEC_WEBSOCKET_KEY1)

        encoding = None
        enc = self._content_encoding
        if enc is not None:
            self._content_encoding = None
            if enc.isascii() and enc.lower() in {"gzip", "deflate", "br", "zstd"}:
                encoding = enc

        if self._cparser.type == cparser.HTTP_REQUEST:
            method = <str>_http_method[self._cparser.method]
            msg = _new_request_message(
                method, self._path,
                http_version, headers, raw_headers,
                should_close, encoding, upgrade, chunked, self._url)
        else:
            msg = _new_response_message(
                http_version, self._cparser.status_code, self._reason,
                headers, raw_headers, should_close, encoding,
                upgrade, chunked)

        if (
            self._response_with_body
            and (
                ULLONG_MAX > self._cparser.content_length > 0 or chunked or
                self._cparser.method == cparser.HTTP_CONNECT or
                (self._cparser.status_code >= 199 and
                 self._cparser.content_length == 0 and
                 self._read_until_eof)
            )
        ):
            payload = StreamReader(
                self.protocol, timer=self._timer, loop=self._loop,
                limit=self._limit)
        else:
            payload = EMPTY_PAYLOAD

        self._payload = payload
        self._content_length_expected = self._cparser.content_length
        if encoding is not None and self._auto_decompress:
            self._payload = DeflateBuffer(payload, encoding, max_decompress_size=self._limit)

        self._messages.append((msg, payload))

    cdef _on_message_complete(self):
        self._payload.feed_eof()
        self._payload = None

    cdef _on_chunk_header(self):
        self._payload.begin_http_chunk_receiving()

    cdef _on_chunk_complete(self):
        self._payload.end_http_chunk_receiving()

    cdef object _on_status_complete(self):
        pass

    cdef inline http_version(self):
        cdef cparser.llhttp_t* parser = self._cparser

        if parser.http_major == 1:
            if parser.http_minor == 0:
                return HttpVersion10
            elif parser.http_minor == 1:
                return HttpVersion11

        return HttpVersion(parser.http_major, parser.http_minor)

    ### Public API ###

    def pause_reading(self):
        assert self._payload is not None
        self._paused = True

    def message_consumed(self):
        # Protocol drained a queued message; free a slot for parsing.
        if self._msg_in_flight > 0:
            self._msg_in_flight -= 1

    def feed_eof(self):
        cdef bytes desc

        if self._payload is not None:
            if self._cparser.flags & cparser.F_CHUNKED:
                raise TransferEncodingError(
                    "Not enough data to satisfy transfer length header.")
            elif self._cparser.flags & cparser.F_CONTENT_LENGTH:
                received = self._content_length_expected - self._cparser.content_length
                raise ContentLengthError(
                    f"Not enough data to satisfy content length header "
                    f"(received {received} of {self._content_length_expected} bytes).")
            elif cparser.llhttp_get_errno(self._cparser) != cparser.HPE_OK:
                desc = cparser.llhttp_get_error_reason(self._cparser)
                raise PayloadEncodingError(desc.decode('latin-1'))
            else:
                self._eof_pending = True
                while self._more_data_available:
                    if self._paused:
                        self._paused = False
                        return  # Will resume via feed_data(b"") later
                    self._more_data_available = self._payload.feed_data(b"", 0)
                self._payload.feed_eof()
                self._payload = None
                self._more_data_available = False
                self._eof_pending = False
        elif self._started:
            self._on_headers_complete()
            if self._messages:
                return self._messages[-1][0]

    def feed_data(self, incoming_data):
        cdef:
            size_t data_len
            size_t nb
            char* base
            cdef cparser.llhttp_errno_t errno
            cdef bytes data

        # Proactor loop sends bytearray.
        # Ensure cython sees `data` as bytes
        if type(incoming_data) is not bytes:
            data = bytes(incoming_data)
        else:
            data = incoming_data

        if self._tail:
            data, self._tail = self._tail + data, b""

        if self._more_data_available:
            result = cb_on_body(self._cparser, b"", 0)
            if result is cparser.HPE_PAUSED:
                self._tail = data
                return EMPTY_FEED_DATA_RESULT

        if self._eof_pending:
            self._payload.feed_eof()
            self._payload = None
            self._eof_pending = False
            # We can't have new messages here, otherwise we wouldn't have
            # received EOF.
            return EMPTY_FEED_DATA_RESULT

        PyObject_GetBuffer(data, &self.py_buf, PyBUF_SIMPLE)
        # Cache buffer pointer before PyBuffer_Release to avoid use-after-release.
        base = <char*>self.py_buf.buf
        data_len = <size_t>self.py_buf.len

        errno = cparser.llhttp_execute(
            self._cparser,
            base,
            data_len)

        if errno is cparser.HPE_PAUSED_UPGRADE:
            cparser.llhttp_resume_after_upgrade(self._cparser)
            nb = cparser.llhttp_get_error_pos(self._cparser) - base
            if self._pending_upgrade:
                # A supported upgrade whose request body has now been fully read.
                self._upgraded = True
                self._pending_upgrade = False
        elif errno is cparser.HPE_PAUSED:
            cparser.llhttp_resume(self._cparser)
            pos = cparser.llhttp_get_error_pos(self._cparser) - base
            self._tail = data[pos:]

        PyBuffer_Release(&self.py_buf)

        if errno not in (cparser.HPE_OK, cparser.HPE_PAUSED, cparser.HPE_PAUSED_UPGRADE):
            if self._payload_error == 0:
                if self._last_error is not None:
                    ex = self._last_error
                    self._last_error = None
                else:
                    after = cparser.llhttp_get_error_pos(self._cparser)
                    before = data[:after - base]
                    after_b = after.split(b"\r\n", 1)[0]
                    before = before.rsplit(b"\r\n", 1)[-1]
                    data = before + after_b
                    pointer = " " * (len(repr(before))-1) + "^"
                    ex = parser_error_from_errno(self._cparser, data, pointer)
                self._payload = None
                raise ex

        if self._messages:
            messages = self._messages
            self._messages = []
        else:
            messages = ()

        if self._upgraded:
            return messages, True, data[nb:]
        if not messages:  # Shortcut to reduce Python overhead
            return EMPTY_FEED_DATA_RESULT
        return messages, False, b""

    def set_upgraded(self, val):
        self._upgraded = val


cdef class HttpRequestParser(HttpParser):

    def __init__(
        self, protocol, loop, int limit, timer=None,
        size_t max_line_size=8190, size_t max_headers=128,
        size_t max_field_size=8190, payload_exception=None,
        bint response_with_body=True, bint read_until_eof=False,
        bint auto_decompress=True, Py_ssize_t max_msg_queue_size=0,
    ):
        self._init(cparser.HTTP_REQUEST, protocol, loop, limit, timer,
                   max_line_size, max_headers, max_field_size,
                   payload_exception, response_with_body, read_until_eof,
                   auto_decompress, max_msg_queue_size)

    cdef object _on_status_complete(self):
        cdef int idx1, idx2
        if not self._buf:
            return
        self._path = self._buf.decode('utf-8', 'surrogateescape')
        try:
            idx3 = len(self._path)
            if self._cparser.method == cparser.HTTP_CONNECT:
                # authority-form,
                # https://datatracker.ietf.org/doc/html/rfc7230#section-5.3.3
                self._url = URL.build(authority=self._path, encoded=True)
            elif idx3 > 1 and self._path[0] == '/':
                # origin-form,
                # https://datatracker.ietf.org/doc/html/rfc7230#section-5.3.1
                idx1 = self._path.find("?")
                if idx1 == -1:
                    query = ""
                    idx2 = self._path.find("#")
                    if idx2 == -1:
                        path = self._path
                        fragment = ""
                    else:
                        path = self._path[0: idx2]
                        fragment = self._path[idx2+1:]

                else:
                    path = self._path[0:idx1]
                    idx1 += 1
                    idx2 = self._path.find("#", idx1)
                    if idx2 == -1:
                        query = self._path[idx1:]
                        fragment = ""
                    else:
                        query = self._path[idx1: idx2]
                        fragment = self._path[idx2+1:]

                self._url = URL.build(
                    path=path,
                    query_string=query,
                    fragment=fragment,
                    encoded=True,
                )
            else:
                # absolute-form for proxy maybe,
                # https://datatracker.ietf.org/doc/html/rfc7230#section-5.3.2
                self._url = URL(self._path, encoded=True)
        finally:
            PyByteArray_Resize(self._buf, 0)


cdef class HttpResponseParser(HttpParser):

    def __init__(
        self, protocol, loop, int limit, timer=None,
            size_t max_line_size=8190, size_t max_headers=128,
            size_t max_field_size=8190, payload_exception=None,
            bint response_with_body=True, bint read_until_eof=False,
            bint auto_decompress=True
    ):
        self._init(cparser.HTTP_RESPONSE, protocol, loop, limit, timer,
                   max_line_size, max_headers, max_field_size,
                   payload_exception, response_with_body, read_until_eof,
                   auto_decompress)
        # Use strict parsing on dev mode, so users are warned about broken servers.
        if not DEBUG:
            cparser.llhttp_set_lenient_headers(self._cparser, 1)
            cparser.llhttp_set_lenient_optional_cr_before_lf(self._cparser, 1)
            cparser.llhttp_set_lenient_spaces_after_chunk_size(self._cparser, 1)
            self._lax = True

    cdef object _on_status_complete(self):
        if self._buf:
            self._reason = self._buf.decode('utf-8', 'surrogateescape')
            PyByteArray_Resize(self._buf, 0)
        else:
            self._reason = self._reason or ''

cdef int cb_on_message_begin(cparser.llhttp_t* parser) except -1:
    cdef HttpParser pyparser = <HttpParser>parser.data

    pyparser._started = True
    pyparser._headers = []
    pyparser._seen_singletons = set()
    pyparser._raw_headers = []
    PyByteArray_Resize(pyparser._buf, 0)
    pyparser._path = None
    pyparser._reason = None
    return 0


cdef int cb_on_url(cparser.llhttp_t* parser,
                   const char *at, size_t length) except -1:
    cdef HttpParser pyparser = <HttpParser>parser.data
    try:
        if len(pyparser._buf) + length > pyparser._max_line_size:
            status = pyparser._buf + at[:length]
            raise LineTooLong(status[:100] + b"...", pyparser._max_line_size)
        extend(pyparser._buf, at, length)
    except BaseException as ex:
        pyparser._last_error = ex
        return -1
    else:
        return 0


cdef int cb_on_status(cparser.llhttp_t* parser,
                      const char *at, size_t length) except -1:
    cdef HttpParser pyparser = <HttpParser>parser.data
    try:
        if len(pyparser._buf) + length > pyparser._max_line_size:
            reason = pyparser._buf + at[:length]
            raise LineTooLong(reason[:100] + b"...", pyparser._max_line_size)
        extend(pyparser._buf, at, length)
    except BaseException as ex:
        pyparser._last_error = ex
        return -1
    else:
        return 0


cdef int cb_on_header_field(cparser.llhttp_t* parser,
                            const char *at, size_t length) except -1:
    cdef HttpParser pyparser = <HttpParser>parser.data
    cdef Py_ssize_t size
    try:
        pyparser._on_status_complete()
        size = len(pyparser._raw_name) + length
        if size > pyparser._max_field_size:
            name = pyparser._raw_name + at[:length]
            raise LineTooLong(name[:100] + b"...", pyparser._max_field_size)
        pyparser._header_name_size = size
        pyparser._on_header_field(at, length)
    except BaseException as ex:
        pyparser._last_error = ex
        return -1
    else:
        return 0


cdef int cb_on_header_value(cparser.llhttp_t* parser,
                            const char *at, size_t length) except -1:
    cdef HttpParser pyparser = <HttpParser>parser.data
    cdef Py_ssize_t size
    try:
        size = len(pyparser._raw_value) + length
        if pyparser._header_name_size + size > pyparser._max_field_size:
            value = pyparser._raw_value + at[:length]
            raise LineTooLong(value[:100] + b"...", pyparser._max_field_size)
        pyparser._on_header_value(at, length)
    except BaseException as ex:
        pyparser._last_error = ex
        return -1
    else:
        return 0


cdef int cb_on_headers_complete(cparser.llhttp_t* parser) except -1:
    cdef HttpParser pyparser = <HttpParser>parser.data
    try:
        pyparser._on_status_complete()
        pyparser._on_headers_complete()
    except BaseException as exc:
        pyparser._last_error = exc
        return -1
    else:
        if not pyparser._response_with_body:
            return 1
        return 0


cdef int cb_on_body(cparser.llhttp_t* parser,
                    const char *at, size_t length) except -1:
    cdef HttpParser pyparser = <HttpParser>parser.data
    cdef bytes body = at[:length]
    while body or pyparser._more_data_available:
        try:
            pyparser._more_data_available = pyparser._payload.feed_data(body, length)
        except BaseException as underlying_exc:
            reraised_exc = underlying_exc
            if pyparser._payload_exception is not None:
                reraised_exc = pyparser._payload_exception(str(underlying_exc))

            set_exception(pyparser._payload, reraised_exc, underlying_exc)

            pyparser._payload_error = 1
            pyparser._paused = False
            return -1
        body = b""
        length = 0

        if pyparser._paused:
            pyparser._paused = False
            return cparser.HPE_PAUSED
    pyparser._paused = False
    return 0


cdef int cb_on_message_complete(cparser.llhttp_t* parser) except -1:
    cdef HttpParser pyparser = <HttpParser>parser.data
    try:
        pyparser._started = False
        pyparser._on_message_complete()
    except BaseException as exc:
        pyparser._last_error = exc
        return -1
    else:
        if pyparser._max_msg_queue_size:
            pyparser._msg_in_flight += 1
            if pyparser._msg_in_flight >= pyparser._max_msg_queue_size:
                # Queue full: pause llhttp between messages. feed_data() buffers
                # the remainder as tail; resumes once the queue drains.
                return cparser.HPE_PAUSED
        return 0


cdef int cb_on_chunk_header(cparser.llhttp_t* parser) except -1:
    cdef HttpParser pyparser = <HttpParser>parser.data
    try:
        pyparser._on_chunk_header()
    except BaseException as exc:
        pyparser._last_error = exc
        return -1
    else:
        return 0


cdef int cb_on_chunk_complete(cparser.llhttp_t* parser) except -1:
    cdef HttpParser pyparser = <HttpParser>parser.data
    try:
        pyparser._on_chunk_complete()
    except BaseException as exc:
        pyparser._last_error = exc
        return -1
    else:
        return 0


cdef parser_error_from_errno(cparser.llhttp_t* parser, data, pointer):
    cdef cparser.llhttp_errno_t errno = cparser.llhttp_get_errno(parser)
    cdef bytes desc = cparser.llhttp_get_error_reason(parser)

    err_msg = "{}:\n\n  {!r}\n  {}".format(desc.decode("latin-1"), data, pointer)

    if errno in {cparser.HPE_CB_MESSAGE_BEGIN,
                 cparser.HPE_CB_HEADERS_COMPLETE,
                 cparser.HPE_CB_MESSAGE_COMPLETE,
                 cparser.HPE_CB_CHUNK_HEADER,
                 cparser.HPE_CB_CHUNK_COMPLETE,
                 cparser.HPE_INVALID_HEADER_TOKEN,
                 cparser.HPE_INVALID_CONTENT_LENGTH,
                 cparser.HPE_INVALID_CHUNK_SIZE,
                 cparser.HPE_INVALID_EOF_STATE,
                 cparser.HPE_INVALID_TRANSFER_ENCODING}:
        return BadHttpMessage(err_msg)
    elif errno == cparser.HPE_INVALID_METHOD:
        if data.startswith(b"\x16\x03"):
            return BadHttpMethod(error="Received HTTPS traffic on an HTTP port")
        return BadHttpMethod(error=err_msg)
    elif errno in {cparser.HPE_INVALID_STATUS,
                   cparser.HPE_INVALID_VERSION,
                   cparser.HPE_INVALID_CONSTANT}:
        return BadStatusLine(error=f"Bad status line:\n  {err_msg}")
    elif errno == cparser.HPE_INVALID_URL:
        return InvalidURLError(err_msg)

    return BadHttpMessage(err_msg)
