import abc
import asyncio
import re
import string
import sys
from contextlib import suppress
from enum import IntEnum
from re import Pattern
from typing import (
    TYPE_CHECKING,
    Any,
    ClassVar,
    Final,
    Generic,
    Literal,
    NamedTuple,
    TypeVar,
)

from multidict import CIMultiDict, CIMultiDictProxy, istr
from yarl import URL

from . import hdrs
from .base_protocol import BaseProtocol
from .compression_utils import (
    HAS_BROTLI,
    HAS_ZSTD,
    BrotliDecompressor,
    ZLibDecompressor,
    ZSTDDecompressor,
)
from .helpers import (
    _EXC_SENTINEL,
    DEBUG,
    DEFAULT_CHUNK_SIZE,
    EMPTY_BODY_METHODS,
    EMPTY_BODY_STATUS_CODES,
    NO_EXTENSIONS,
    BaseTimerContext,
    set_exception,
)
from .http_exceptions import (
    BadHttpMessage,
    BadHttpMethod,
    BadStatusLine,
    ContentEncodingError,
    ContentLengthError,
    InvalidHeader,
    InvalidURLError,
    LineTooLong,
    TransferEncodingError,
)
from .http_writer import HttpVersion, HttpVersion10, HttpVersion11
from .streams import EMPTY_PAYLOAD, StreamReader
from .typedefs import RawHeaders

if TYPE_CHECKING:
    from .client_proto import ResponseHandler

__all__ = (
    "HeadersParser",
    "HttpParser",
    "HttpRequestParser",
    "HttpResponseParser",
    "RawRequestMessage",
    "RawResponseMessage",
)

_SEP = Literal[b"\r\n", b"\n"]

ASCIISET: Final[set[str]] = set(string.printable)

# See https://www.rfc-editor.org/rfc/rfc9110.html#name-overview
# and https://www.rfc-editor.org/rfc/rfc9110.html#name-tokens
#
#     method = token
#     tchar = "!" / "#" / "$" / "%" / "&" / "'" / "*" / "+" / "-" / "." /
#             "^" / "_" / "`" / "|" / "~" / DIGIT / ALPHA
#     token = 1*tchar
_TCHAR_SPECIALS: Final[str] = re.escape("!#$%&'*+-.^_`|~")
TOKENRE: Final[Pattern[str]] = re.compile(f"[0-9A-Za-z{_TCHAR_SPECIALS}]+")
VERSRE: Final[Pattern[str]] = re.compile(r"HTTP/(\d)\.(\d)", re.ASCII)
DIGITS: Final[Pattern[str]] = re.compile(r"\d+", re.ASCII)
HEXDIGITS: Final[Pattern[bytes]] = re.compile(rb"[0-9a-fA-F]+")
# https://www.rfc-editor.org/rfc/rfc9110#section-5.5-5
_FIELD_VALUE_FORBIDDEN_CTL_RE: Final[Pattern[str]] = re.compile(
    r"[\x00-\x08\x0a-\x1f\x7f]"
)

# RFC 9110 singleton headers — duplicates are rejected in strict mode.
# In lax mode (response parser default), the check is skipped entirely
# since real-world servers (e.g. Google APIs, Werkzeug) commonly send
# duplicate headers like Content-Type or Server.
# Lowercased for case-insensitive matching against wire names.
SINGLETON_HEADERS: Final[frozenset[str]] = frozenset(
    {
        "content-length",
        "content-location",
        "content-range",
        "content-type",
        "etag",
        "host",
        "max-forwards",
        "server",
        "transfer-encoding",
        "user-agent",
    }
)


class RawRequestMessage(NamedTuple):
    method: str
    path: str
    version: HttpVersion
    headers: "CIMultiDictProxy[str]"
    raw_headers: RawHeaders
    should_close: bool
    compression: str | None
    upgrade: bool
    chunked: bool
    url: URL


class RawResponseMessage(NamedTuple):
    version: HttpVersion
    code: int
    reason: str
    headers: CIMultiDictProxy[str]
    raw_headers: RawHeaders
    should_close: bool
    compression: str | None
    upgrade: bool
    chunked: bool


_MsgT = TypeVar("_MsgT", RawRequestMessage, RawResponseMessage)


class PayloadState(IntEnum):
    PAYLOAD_COMPLETE = 0
    PAYLOAD_NEEDS_INPUT = 1
    PAYLOAD_HAS_PENDING_INPUT = 2


class ParseState(IntEnum):

    PARSE_NONE = 0
    PARSE_LENGTH = 1
    PARSE_CHUNKED = 2
    PARSE_UNTIL_EOF = 3


class ChunkState(IntEnum):
    PARSE_CHUNKED_SIZE = 0
    PARSE_CHUNKED_CHUNK = 1
    PARSE_CHUNKED_CHUNK_EOF = 2
    PARSE_MAYBE_TRAILERS = 3
    PARSE_TRAILERS = 4


class HeadersParser:
    def __init__(
        self,
        max_line_size: int = 8190,
        max_headers: int = 32768,
        max_field_size: int = 8190,
        lax: bool = False,
    ) -> None:
        self.max_line_size = max_line_size
        self.max_headers = max_headers
        self.max_field_size = max_field_size
        self._lax = lax

    def parse_headers(
        self, lines: list[bytes]
    ) -> tuple["CIMultiDictProxy[str]", RawHeaders]:
        headers: CIMultiDict[str] = CIMultiDict()
        # note: "raw" does not mean inclusion of OWS before/after the field value
        raw_headers = []

        lines_idx = 0
        line = lines[lines_idx]
        line_count = len(lines)

        while line:
            # Parse initial header name : value pair.
            try:
                bname, bvalue = line.split(b":", 1)
            except ValueError:
                raise InvalidHeader(line) from None

            if len(bname) == 0:
                raise InvalidHeader(bname)

            # https://www.rfc-editor.org/rfc/rfc9112.html#section-5.1-2
            if {bname[0], bname[-1]} & {32, 9}:  # {" ", "\t"}
                raise InvalidHeader(line)

            bvalue = bvalue.lstrip(b" \t")
            name = bname.decode("utf-8", "surrogateescape")
            if not TOKENRE.fullmatch(name):
                raise InvalidHeader(bname)

            # next line
            lines_idx += 1
            line = lines[lines_idx]

            # consume continuation lines
            continuation = self._lax and line and line[0] in (32, 9)  # (' ', '\t')

            # Deprecated: https://www.rfc-editor.org/rfc/rfc9112.html#name-obsolete-line-folding
            if continuation:
                header_length = len(bvalue)
                bvalue_lst = [bvalue]
                while continuation:
                    header_length += len(line)
                    if header_length > self.max_field_size:
                        header_line = bname + b": " + b"".join(bvalue_lst)
                        raise LineTooLong(
                            header_line[:100] + b"...", self.max_field_size
                        )
                    bvalue_lst.append(line)

                    # next line
                    lines_idx += 1
                    if lines_idx < line_count:
                        line = lines[lines_idx]
                        if line:
                            continuation = line[0] in (32, 9)  # (' ', '\t')
                    else:
                        line = b""
                        break
                bvalue = b"".join(bvalue_lst)

            bvalue = bvalue.strip(b" \t")
            value = bvalue.decode("utf-8", "surrogateescape")

            # https://www.rfc-editor.org/rfc/rfc9110.html#section-5.5-5
            if self._lax:
                if "\n" in value or "\r" in value or "\x00" in value:
                    raise InvalidHeader(bvalue)
            elif _FIELD_VALUE_FORBIDDEN_CTL_RE.search(value):
                raise InvalidHeader(bvalue)

            if not self._lax and name in headers and name.lower() in SINGLETON_HEADERS:
                raise BadHttpMessage(f"Duplicate '{name}' header found.")
            headers.add(name, value)
            raw_headers.append((bname, bvalue))

        return (CIMultiDictProxy(headers), tuple(raw_headers))


def _is_supported_upgrade(headers: CIMultiDictProxy[str]) -> bool:
    """Check if the upgrade header is supported."""
    u = headers.get(hdrs.UPGRADE, "")
    # .lower() can transform non-ascii characters.
    return u.isascii() and u.lower() in {"tcp", "websocket"}


class HttpParser(abc.ABC, Generic[_MsgT]):
    lax: ClassVar[bool] = False

    def __init__(
        self,
        protocol: BaseProtocol | None = None,
        loop: asyncio.AbstractEventLoop | None = None,
        limit: int = 2**16,
        max_line_size: int = 8190,
        max_headers: int = 128,
        max_field_size: int = 8190,
        timer: BaseTimerContext | None = None,
        code: int | None = None,
        method: str | None = None,
        payload_exception: type[BaseException] | None = None,
        response_with_body: bool = True,
        read_until_eof: bool = False,
        auto_decompress: bool = True,
        max_msg_queue_size: int = 0,
    ) -> None:
        self.protocol = protocol
        self.loop = loop
        self.max_line_size = max_line_size
        self.max_headers = max_headers
        self.max_field_size = max_field_size
        self.max_headers = max_headers
        self.timer = timer
        self.code = code
        self.method = method
        self.payload_exception = payload_exception
        self.response_with_body = response_with_body
        self.read_until_eof = read_until_eof

        self._lines: list[bytes] = []
        self._tail = b""
        self._upgraded = False
        self._pending_upgrade = False
        self._payload = None
        self._payload_parser: HttpPayloadParser | None = None
        self._payload_has_more_data = False
        self._auto_decompress = auto_decompress
        self._limit = limit
        self._headers_parser = HeadersParser(
            max_line_size, max_headers, max_field_size, self.lax
        )
        # Stop emitting messages once this many are queued unconsumed (0 = off).
        self._max_msg_queue_size = max_msg_queue_size
        self._msg_in_flight = 0

    @abc.abstractmethod
    def parse_message(self, lines: list[bytes]) -> _MsgT: ...

    @abc.abstractmethod
    def _is_chunked_te(self, te: str) -> bool: ...

    def pause_reading(self) -> None:
        assert self._payload_parser is not None
        self._payload_parser.pause_reading()

    def message_consumed(self) -> None:
        """Protocol drained a queued message; free a slot for parsing."""
        if self._msg_in_flight > 0:
            self._msg_in_flight -= 1

    def feed_eof(self) -> _MsgT | None:
        if self._payload_parser is not None:
            self._payload_parser.feed_eof()
            if self._payload_parser.done:
                self._payload_parser = None
        else:
            # try to extract partial message
            if self._tail:
                self._lines.append(self._tail)

            if self._lines:
                if self._lines[-1] != "\r\n":
                    self._lines.append(b"")
                with suppress(Exception):
                    return self.parse_message(self._lines)
        return None

    def feed_data(
        self,
        data: bytes,
        SEP: _SEP = b"\r\n",
        EMPTY: bytes = b"",
        CONTENT_LENGTH: istr = hdrs.CONTENT_LENGTH,
        METH_CONNECT: str = hdrs.METH_CONNECT,
        SEC_WEBSOCKET_KEY1: istr = hdrs.SEC_WEBSOCKET_KEY1,
    ) -> tuple[list[tuple[_MsgT, StreamReader]], bool, bytes]:

        messages = []

        if self._tail:
            data, self._tail = self._tail + data, b""

        data_len = len(data)
        start_pos = 0
        loop = self.loop
        max_line_length = self.max_line_size

        should_close = False
        while start_pos < data_len or self._payload_has_more_data:
            # read HTTP message (request/response line + headers), \r\n\r\n
            # and split by lines
            if self._payload_parser is None and not self._upgraded:
                if (
                    self._max_msg_queue_size
                    and self._msg_in_flight >= self._max_msg_queue_size
                ):
                    # Queue full: buffer the rest and stop. Safe pause point;
                    # any preceding body is consumed before the next request
                    # line. Resumes via feed_data(b"") when the queue drains.
                    self._tail = data[start_pos:]
                    break
                pos = data.find(SEP, start_pos)
                # consume \r\n
                if pos == start_pos and not self._lines:
                    start_pos = pos + len(SEP)
                    continue

                if pos >= start_pos:
                    if should_close:
                        raise BadHttpMessage("Data after `Connection: close`")

                    # line found
                    line = data[start_pos:pos]
                    if SEP == b"\n":  # For lax response parsing
                        line = line.rstrip(b"\r")
                    if len(line) > max_line_length:
                        raise LineTooLong(line[:100] + b"...", max_line_length)

                    self._lines.append(line)
                    # After processing the status/request line, everything is a header.
                    max_line_length = self.max_field_size

                    if len(self._lines) > self.max_headers:
                        raise BadHttpMessage("Too many headers received")

                    start_pos = pos + len(SEP)

                    # \r\n\r\n found
                    if self._lines[-1] == EMPTY:
                        max_trailers = self.max_headers - len(self._lines)
                        try:
                            msg: _MsgT = self.parse_message(self._lines)
                        finally:
                            self._lines.clear()

                        def get_content_length() -> int | None:
                            # payload length
                            length_hdr = msg.headers.get(CONTENT_LENGTH)
                            if length_hdr is None:
                                return None

                            # Shouldn't allow +/- or other number formats.
                            # https://www.rfc-editor.org/rfc/rfc9110#section-8.6-2
                            # msg.headers is already stripped of leading/trailing wsp
                            if not DIGITS.fullmatch(length_hdr):
                                raise InvalidHeader(CONTENT_LENGTH)

                            return int(length_hdr)

                        length = get_content_length()
                        # do not support old websocket spec
                        if SEC_WEBSOCKET_KEY1 in msg.headers:
                            raise InvalidHeader(SEC_WEBSOCKET_KEY1)

                        upgraded = msg.upgrade and _is_supported_upgrade(msg.headers)

                        method = getattr(msg, "method", self.method)
                        # code is only present on responses
                        code = getattr(msg, "code", 0)

                        assert self.protocol is not None
                        # calculate payload
                        empty_body = code in EMPTY_BODY_STATUS_CODES or bool(
                            method and method in EMPTY_BODY_METHODS
                        )
                        if not empty_body and (
                            (length is not None and length > 0) or msg.chunked
                        ):
                            payload = StreamReader(
                                self.protocol,
                                timer=self.timer,
                                loop=loop,
                                limit=self._limit,
                            )
                            payload_parser = HttpPayloadParser(
                                payload,
                                length=length,
                                chunked=msg.chunked,
                                method=method,
                                compression=msg.compression,
                                code=self.code,
                                response_with_body=self.response_with_body,
                                auto_decompress=self._auto_decompress,
                                lax=self.lax,
                                headers_parser=self._headers_parser,
                                max_line_size=self.max_line_size,
                                max_field_size=self.max_field_size,
                                max_trailers=max_trailers,
                                limit=self._limit,
                            )
                            if not payload_parser.done:
                                self._payload_parser = payload_parser
                                # https://www.rfc-editor.org/info/rfc9110/#section-7.8-15
                                # Defer any requested upgrade until the
                                # complete request has been read.
                                self._pending_upgrade = upgraded
                        elif method == METH_CONNECT:
                            assert isinstance(msg, RawRequestMessage)
                            payload = StreamReader(
                                self.protocol,
                                timer=self.timer,
                                loop=loop,
                                limit=self._limit,
                            )
                            self._upgraded = True
                            self._payload_parser = HttpPayloadParser(
                                payload,
                                method=msg.method,
                                compression=msg.compression,
                                auto_decompress=self._auto_decompress,
                                lax=self.lax,
                                headers_parser=self._headers_parser,
                                max_line_size=self.max_line_size,
                                max_field_size=self.max_field_size,
                                max_trailers=max_trailers,
                                limit=self._limit,
                            )
                        elif not empty_body and length is None and self.read_until_eof:
                            payload = StreamReader(
                                self.protocol,
                                timer=self.timer,
                                loop=loop,
                                limit=self._limit,
                            )
                            payload_parser = HttpPayloadParser(
                                payload,
                                length=length,
                                chunked=msg.chunked,
                                method=method,
                                compression=msg.compression,
                                code=self.code,
                                response_with_body=self.response_with_body,
                                auto_decompress=self._auto_decompress,
                                lax=self.lax,
                                headers_parser=self._headers_parser,
                                max_line_size=self.max_line_size,
                                max_field_size=self.max_field_size,
                                max_trailers=max_trailers,
                                limit=self._limit,
                            )
                            if not payload_parser.done:
                                self._payload_parser = payload_parser
                        elif upgraded:
                            # No body to read, so the connection switches to
                            # the upgraded protocol immediately.
                            self._upgraded = True
                            payload = EMPTY_PAYLOAD
                        else:
                            payload = EMPTY_PAYLOAD

                        messages.append((msg, payload))
                        if self._max_msg_queue_size:
                            self._msg_in_flight += 1
                        should_close = msg.should_close
                else:
                    self._tail = data[start_pos:]
                    # A bare LF here means CRLF was required:
                    # reject instead of buffering, else a following request's
                    # bytes get appended to this line and leak in the error.
                    if b"\n" in self._tail:
                        raise BadHttpMessage("Bad line ending, expected CRLF")
                    if len(self._tail) > self.max_line_size:
                        raise LineTooLong(self._tail[:100] + b"...", self.max_line_size)
                    data = EMPTY
                    break

            # no parser, just store
            elif self._payload_parser is None and self._upgraded:
                assert not self._lines
                break

            # feed payload
            else:
                assert not self._lines
                assert self._payload_parser is not None
                try:
                    payload_state, data = self._payload_parser.feed_data(
                        data[start_pos:], SEP
                    )
                except Exception as underlying_exc:
                    reraised_exc: BaseException = underlying_exc
                    if self.payload_exception is not None:
                        reraised_exc = self.payload_exception(str(underlying_exc))

                    set_exception(
                        self._payload_parser.payload,
                        reraised_exc,
                        underlying_exc,
                    )

                    payload_state = PayloadState.PAYLOAD_COMPLETE
                    data = b""
                    if isinstance(
                        underlying_exc, (InvalidHeader, TransferEncodingError)
                    ):
                        raise

                self._payload_has_more_data = (
                    payload_state == PayloadState.PAYLOAD_HAS_PENDING_INPUT
                )

                if payload_state is not PayloadState.PAYLOAD_COMPLETE:
                    # We've either consumed all available data, or we're pausing
                    # until the reader buffer is freed up.
                    break

                start_pos = 0
                data_len = len(data)
                self._payload_parser = None
                if self._pending_upgrade:
                    # Body fully read: the deferred upgrade takes effect and
                    # the rest of the connection is the upgraded protocol.
                    self._upgraded = True
                    self._pending_upgrade = False

        if data and start_pos < data_len:
            data = data[start_pos:]
        else:
            data = EMPTY

        return messages, self._upgraded, data

    def parse_headers(
        self, lines: list[bytes]
    ) -> tuple[
        "CIMultiDictProxy[str]", RawHeaders, bool | None, str | None, bool, bool
    ]:
        """Parses RFC 5322 headers from a stream.

        Line continuations are supported. Returns list of header name
        and value pairs. Header name is in upper case.
        """
        headers, raw_headers = self._headers_parser.parse_headers(lines)
        close_conn = None
        encoding = None
        upgrade = False
        chunked = False

        # keep-alive and protocol switching
        # RFC 9110 section 7.6.1 defines Connection as a comma-separated list.
        conn_values = headers.getall(hdrs.CONNECTION, ())
        if conn_values:
            conn_tokens = {
                token.lower()
                for conn_value in conn_values
                for token in (part.strip(" \t") for part in conn_value.split(","))
                if token and token.isascii()
            }

            if "close" in conn_tokens:
                close_conn = True
            elif "keep-alive" in conn_tokens:
                close_conn = False

            # https://www.rfc-editor.org/rfc/rfc9110.html#name-101-switching-protocols
            if "upgrade" in conn_tokens and headers.get(hdrs.UPGRADE):
                upgrade = True

        # encoding
        enc = headers.get(hdrs.CONTENT_ENCODING, "")
        if enc.isascii() and enc.lower() in {"gzip", "deflate", "br", "zstd"}:
            encoding = enc

        # chunking
        te = headers.get(hdrs.TRANSFER_ENCODING)
        if te is not None:
            if self._is_chunked_te(te):
                chunked = True

            if hdrs.CONTENT_LENGTH in headers:
                raise BadHttpMessage(
                    "Transfer-Encoding can't be present with Content-Length",
                )

        return (headers, raw_headers, close_conn, encoding, upgrade, chunked)

    def set_upgraded(self, val: bool) -> None:
        """Set connection upgraded (to websocket) mode.

        :param bool val: new state.
        """
        self._upgraded = val


class HttpRequestParser(HttpParser[RawRequestMessage]):
    """Read request status line.

    Exception .http_exceptions.BadStatusLine
    could be raised in case of any errors in status line.
    Returns RawRequestMessage.
    """

    def parse_message(self, lines: list[bytes]) -> RawRequestMessage:
        # request line
        line = lines[0].decode("utf-8", "surrogateescape")
        try:
            method, path, version = line.split(" ", maxsplit=2)
        except ValueError:
            raise BadHttpMethod(line) from None

        # method
        if not TOKENRE.fullmatch(method):
            raise BadHttpMethod(method)
        method = method.upper()

        # version
        match = VERSRE.fullmatch(version)
        if match is None:
            raise BadStatusLine(line)
        version_o = HttpVersion(int(match.group(1)), int(match.group(2)))

        if method == "CONNECT":
            # authority-form,
            # https://datatracker.ietf.org/doc/html/rfc7230#section-5.3.3
            url = URL.build(authority=path, encoded=True)
        elif path.startswith("/"):
            # origin-form,
            # https://datatracker.ietf.org/doc/html/rfc7230#section-5.3.1
            path_part, _hash_separator, url_fragment = path.partition("#")
            path_part, _question_mark_separator, qs_part = path_part.partition("?")

            # NOTE: `yarl.URL.build()` is used to mimic what the Cython-based
            # NOTE: parser does, otherwise it results into the same
            # NOTE: HTTP Request-Line input producing different
            # NOTE: `yarl.URL()` objects
            url = URL.build(
                path=path_part,
                query_string=qs_part,
                fragment=url_fragment,
                encoded=True,
            )
        elif path == "*" and method == "OPTIONS":
            # asterisk-form,
            url = URL(path, encoded=True)
        else:
            # absolute-form for proxy maybe,
            # https://datatracker.ietf.org/doc/html/rfc7230#section-5.3.2
            url = URL(path, encoded=True)
            if url.scheme == "":
                # not absolute-form
                raise InvalidURLError(
                    path.encode(errors="surrogateescape").decode("latin1")
                )

        # read headers
        (
            headers,
            raw_headers,
            close,
            compression,
            upgrade,
            chunked,
        ) = self.parse_headers(lines[1:])

        if version_o == HttpVersion11 and hdrs.HOST not in headers:
            raise BadHttpMessage("Missing 'Host' header in request.")

        if close is None:  # then the headers weren't set in the request
            if version_o <= HttpVersion10:  # HTTP 1.0 must asks to not close
                close = True
            else:  # HTTP 1.1 must ask to close.
                close = False

        return RawRequestMessage(
            method,
            path,
            version_o,
            headers,
            raw_headers,
            close,
            compression,
            upgrade,
            chunked,
            url,
        )

    def _is_chunked_te(self, te: str) -> bool:
        te = te.rsplit(",", maxsplit=1)[-1].strip(" \t")
        # .lower() transforms some non-ascii chars, so must check first.
        if te.isascii() and te.lower() == "chunked":
            return True
        # https://www.rfc-editor.org/rfc/rfc9112#section-6.3-2.4.3
        raise BadHttpMessage("Request has invalid `Transfer-Encoding`")


class HttpResponseParser(HttpParser[RawResponseMessage]):
    """Read response status line and headers.

    BadStatusLine could be raised in case of any errors in status line.
    Returns RawResponseMessage.
    """

    protocol: "ResponseHandler"

    # Lax mode should only be enabled on response parser.
    lax = not DEBUG

    def feed_data(
        self,
        data: bytes,
        SEP: _SEP | None = None,
        *args: Any,
        **kwargs: Any,
    ) -> tuple[list[tuple[RawResponseMessage, StreamReader]], bool, bytes]:
        if SEP is None:
            SEP = b"\r\n" if DEBUG else b"\n"
        return super().feed_data(data, SEP, *args, **kwargs)

    def parse_message(self, lines: list[bytes]) -> RawResponseMessage:
        line = lines[0].decode("utf-8", "surrogateescape")
        try:
            version, status = line.split(maxsplit=1)
        except ValueError:
            raise BadStatusLine(line) from None

        try:
            status, reason = status.split(maxsplit=1)
        except ValueError:
            status = status.strip()
            reason = ""

        # version
        match = VERSRE.fullmatch(version)
        if match is None:
            raise BadStatusLine(line)
        version_o = HttpVersion(int(match.group(1)), int(match.group(2)))

        # The status code is a three-digit ASCII number, no padding
        if len(status) != 3 or not DIGITS.fullmatch(status):
            raise BadStatusLine(line)
        status_i = int(status)

        # read headers
        (
            headers,
            raw_headers,
            close,
            compression,
            upgrade,
            chunked,
        ) = self.parse_headers(lines[1:])

        if close is None:
            if version_o <= HttpVersion10:
                close = True
            # https://www.rfc-editor.org/rfc/rfc9112.html#name-message-body-length
            elif 100 <= status_i < 200 or status_i in {204, 304}:
                close = False
            elif hdrs.CONTENT_LENGTH in headers or hdrs.TRANSFER_ENCODING in headers:
                close = False
            else:
                # https://www.rfc-editor.org/rfc/rfc9112.html#section-6.3-2.8
                close = True

        return RawResponseMessage(
            version_o,
            status_i,
            reason.strip(),
            headers,
            raw_headers,
            close,
            compression,
            upgrade,
            chunked,
        )

    def _is_chunked_te(self, te: str) -> bool:
        # https://www.rfc-editor.org/rfc/rfc9112#section-6.3-2.4.2
        return te.rsplit(",", maxsplit=1)[-1].strip(" \t").lower() == "chunked"


class HttpPayloadParser:
    def __init__(
        self,
        payload: StreamReader,
        length: int | None = None,
        chunked: bool = False,
        compression: str | None = None,
        code: int | None = None,
        method: str | None = None,
        response_with_body: bool = True,
        auto_decompress: bool = True,
        lax: bool = False,
        *,
        headers_parser: HeadersParser,
        max_line_size: int = 8190,
        max_field_size: int = 8190,
        max_trailers: int = 128,
        limit: int = DEFAULT_CHUNK_SIZE,
    ) -> None:
        self._length = 0
        self._paused = False
        self._type = ParseState.PARSE_UNTIL_EOF
        self._chunk = ChunkState.PARSE_CHUNKED_SIZE
        self._chunk_size = 0
        self._chunk_tail = b""
        self._auto_decompress = auto_decompress
        self._lax = lax
        self._headers_parser = headers_parser
        self._max_line_size = max_line_size
        self._max_field_size = max_field_size
        self._max_trailers = max_trailers
        self._more_data_available = False
        self._trailer_lines: list[bytes] = []
        self.done = False
        self._eof_pending = False

        # payload decompression wrapper
        if response_with_body and compression and self._auto_decompress:
            real_payload: StreamReader | DeflateBuffer = DeflateBuffer(
                payload, compression, max_decompress_size=limit
            )
        else:
            real_payload = payload

        # payload parser
        if not response_with_body:
            # don't parse payload if it's not expected to be received
            self._type = ParseState.PARSE_NONE
            real_payload.feed_eof()
            self.done = True
        elif chunked:
            self._type = ParseState.PARSE_CHUNKED
        elif length is not None:
            self._type = ParseState.PARSE_LENGTH
            self._length = length
            self._length_expected = length
            if self._length == 0:
                real_payload.feed_eof()
                self.done = True

        self.payload = real_payload

    def pause_reading(self) -> None:
        self._paused = True

    def feed_eof(self) -> None:
        if self._type == ParseState.PARSE_UNTIL_EOF:
            self._eof_pending = True
            while self._more_data_available:
                if self._paused:
                    self._paused = False
                    return  # Will resume via feed_data(b"") later
                self._more_data_available = self.payload.feed_data(b"", 0)
            self.payload.feed_eof()
            self.done = True
            self._eof_pending = False
        elif self._type == ParseState.PARSE_LENGTH:
            received = self._length_expected - self._length
            raise ContentLengthError(
                f"Not enough data to satisfy content length header "
                f"(received {received} of {self._length_expected} bytes)."
            )
        elif self._type == ParseState.PARSE_CHUNKED:
            raise TransferEncodingError(
                "Not enough data to satisfy transfer length header."
            )

    def feed_data(
        self, chunk: bytes, SEP: _SEP = b"\r\n", CHUNK_EXT: bytes = b";"
    ) -> tuple[PayloadState, bytes]:
        """Receive a chunk of data to process.

        Return:
            PayloadState - The current state of payload processing.
                           This function may be called with empty bytes after returning
                           PAYLOAD_HAS_PENDING_INPUT to continue processing after a pause.
            bytes - If payload is complete, this is the unconsumed bytes intended for the
                    next message/payload, b"" otherwise.
        """
        # Read specified amount of bytes
        if self._type == ParseState.PARSE_LENGTH:
            if self._chunk_tail:
                chunk = self._chunk_tail + chunk
                self._chunk_tail = b""

            required = self._length
            self._length = max(required - len(chunk), 0)
            self._more_data_available = self.payload.feed_data(
                chunk[:required], required
            )
            while self._more_data_available:
                if self._paused:
                    self._paused = False
                    self._chunk_tail = chunk[required:]
                    return PayloadState.PAYLOAD_HAS_PENDING_INPUT, b""
                self._more_data_available = self.payload.feed_data(b"", 0)

            if self._length == 0:
                self.payload.feed_eof()
                return PayloadState.PAYLOAD_COMPLETE, chunk[required:]
        # Chunked transfer encoding parser
        elif self._type == ParseState.PARSE_CHUNKED:
            if self._chunk_tail:
                # We should check the length is sane when not processing payload body.
                if self._chunk != ChunkState.PARSE_CHUNKED_CHUNK:
                    max_line_length = self._max_line_size
                    if self._chunk == ChunkState.PARSE_TRAILERS:
                        max_line_length = self._max_field_size
                    if len(self._chunk_tail) > max_line_length:
                        raise LineTooLong(
                            self._chunk_tail[:100] + b"...", max_line_length
                        )

                chunk = self._chunk_tail + chunk
                self._chunk_tail = b""

            while chunk or self._more_data_available:
                # read next chunk size
                if self._chunk == ChunkState.PARSE_CHUNKED_SIZE:
                    pos = chunk.find(SEP)
                    if pos >= 0:
                        # Only chunk-size lines reach here; trailers enforce
                        # _max_field_size separately in PARSE_TRAILERS below.
                        if pos > self._max_line_size:
                            raise LineTooLong(chunk[:100] + b"...", self._max_line_size)
                        i = chunk.find(CHUNK_EXT, 0, pos)
                        if i >= 0:
                            size_b = chunk[:i]  # strip chunk-extensions
                            # Verify no LF in the chunk-extension
                            if b"\n" in (ext := chunk[i:pos]):
                                exc = TransferEncodingError(
                                    f"Unexpected LF in chunk-extension: {ext!r}"
                                )
                                set_exception(self.payload, exc)
                                raise exc
                        else:
                            size_b = chunk[:pos]

                        if self._lax:  # Allow whitespace in lax mode.
                            size_b = size_b.strip()

                        if not re.fullmatch(HEXDIGITS, size_b):
                            exc = TransferEncodingError(
                                chunk[:pos].decode("ascii", "surrogateescape")
                            )
                            set_exception(self.payload, exc)
                            raise exc
                        size = int(bytes(size_b), 16)

                        chunk = chunk[pos + len(SEP) :]
                        if size == 0:  # eof marker
                            self._chunk = ChunkState.PARSE_TRAILERS
                            if self._lax and chunk.startswith(b"\r"):
                                chunk = chunk[1:]
                        else:
                            self._chunk = ChunkState.PARSE_CHUNKED_CHUNK
                            self._chunk_size = size
                            self.payload.begin_http_chunk_receiving()
                    else:
                        if b"\n" in chunk:
                            exc = TransferEncodingError(
                                "Bad chunk-size line ending, expected CRLF"
                            )
                            set_exception(self.payload, exc)
                            raise exc
                        self._chunk_tail = chunk
                        return PayloadState.PAYLOAD_NEEDS_INPUT, b""

                # read chunk and feed buffer
                if self._chunk == ChunkState.PARSE_CHUNKED_CHUNK:
                    if self._paused:
                        self._paused = False
                        self._chunk_tail = chunk
                        return PayloadState.PAYLOAD_HAS_PENDING_INPUT, b""

                    required = self._chunk_size
                    self._chunk_size = max(required - len(chunk), 0)
                    self._more_data_available = self.payload.feed_data(
                        chunk[:required], required
                    )
                    chunk = chunk[required:]

                    if self._more_data_available:
                        continue

                    if self._chunk_size:
                        self._paused = False
                        return PayloadState.PAYLOAD_NEEDS_INPUT, b""
                    self._chunk = ChunkState.PARSE_CHUNKED_CHUNK_EOF
                    self.payload.end_http_chunk_receiving()

                # toss the CRLF at the end of the chunk
                if self._chunk == ChunkState.PARSE_CHUNKED_CHUNK_EOF:
                    if self._lax and chunk.startswith(b"\r"):
                        chunk = chunk[1:]
                    if chunk[: len(SEP)] == SEP:
                        chunk = chunk[len(SEP) :]
                        self._chunk = ChunkState.PARSE_CHUNKED_SIZE
                    elif len(chunk) >= len(SEP) or chunk != SEP[: len(chunk)]:
                        exc = TransferEncodingError(
                            "Chunk size mismatch: expected CRLF after chunk data"
                        )
                        set_exception(self.payload, exc)
                        raise exc
                    else:
                        self._chunk_tail = chunk
                        return PayloadState.PAYLOAD_NEEDS_INPUT, b""

                if self._chunk == ChunkState.PARSE_TRAILERS:
                    pos = chunk.find(SEP)
                    if pos < 0:  # No line found
                        if b"\n" in chunk:
                            exc = TransferEncodingError(
                                "Bad trailer line ending, expected CRLF"
                            )
                            set_exception(self.payload, exc)
                            raise exc
                        self._chunk_tail = chunk
                        return PayloadState.PAYLOAD_NEEDS_INPUT, b""

                    line = chunk[:pos]
                    chunk = chunk[pos + len(SEP) :]
                    if SEP == b"\n":  # For lax response parsing
                        line = line.rstrip(b"\r")

                    if len(line) > self._max_field_size:
                        raise LineTooLong(line[:100] + b"...", self._max_field_size)

                    self._trailer_lines.append(line)

                    if len(self._trailer_lines) > self._max_trailers:
                        raise BadHttpMessage("Too many trailers received")

                    # \r\n\r\n found, end of stream
                    if self._trailer_lines[-1] == b"":
                        # Headers and trailers are defined the same way,
                        # so we reuse the HeadersParser here.
                        try:
                            trailers, raw_trailers = self._headers_parser.parse_headers(
                                self._trailer_lines
                            )
                        finally:
                            self._trailer_lines.clear()
                        self.payload.feed_eof()
                        return PayloadState.PAYLOAD_COMPLETE, chunk

        # Read all bytes until eof
        elif self._type == ParseState.PARSE_UNTIL_EOF:
            self._more_data_available = self.payload.feed_data(chunk, len(chunk))
            while self._more_data_available:
                if self._paused:
                    self._paused = False
                    return PayloadState.PAYLOAD_HAS_PENDING_INPUT, b""
                self._more_data_available = self.payload.feed_data(b"", 0)

            if self._eof_pending:
                self.payload.feed_eof()
                self.done = True
                self._eof_pending = False
                return PayloadState.PAYLOAD_COMPLETE, b""

        return PayloadState.PAYLOAD_NEEDS_INPUT, b""


class DeflateBuffer:
    """DeflateStream decompress stream and feed data into specified stream."""

    decompressor: Any

    def __init__(
        self,
        out: StreamReader,
        encoding: str | None,
        max_decompress_size: int = DEFAULT_CHUNK_SIZE,
    ) -> None:
        self.out = out
        self.size = 0
        out.total_compressed_bytes = self.size
        self.encoding = encoding
        self._started_decoding = False

        self.decompressor: BrotliDecompressor | ZLibDecompressor | ZSTDDecompressor
        if encoding == "br":
            if not HAS_BROTLI:  # pragma: no cover
                raise ContentEncodingError(
                    "Can not decode content-encoding: brotli (br). "
                    "Please install `Brotli`"
                )
            self.decompressor = BrotliDecompressor()
        elif encoding == "zstd":
            if not HAS_ZSTD:
                raise ContentEncodingError(
                    "Can not decode content-encoding: zstandard (zstd). "
                    "Please install `backports.zstd`"
                )
            self.decompressor = ZSTDDecompressor()
        else:
            self.decompressor = ZLibDecompressor(encoding=encoding)

        self._max_decompress_size = max_decompress_size

    def set_exception(
        self,
        exc: BaseException,
        exc_cause: BaseException = _EXC_SENTINEL,
    ) -> None:
        set_exception(self.out, exc, exc_cause)

    def feed_data(self, chunk: bytes, size: int) -> bool:
        self.size += size
        self.out.total_compressed_bytes = self.size

        # Inspect the first real byte once to choose the decompressor. An empty
        # chunk (e.g. a chunk-size line arriving without body bytes) has no
        # header to sniff, so skip it and wait for the first data byte.
        if not self._started_decoding and chunk:
            # RFC1950
            # bits 0..3 = CM = 0b1000 = 8 = "deflate"
            # bits 4..7 = CINFO = 1..7 = windows size.
            if self.encoding == "deflate" and chunk[0] & 0xF != 8:
                # Change the decoder to decompress incorrectly compressed data
                # Actually we should issue a warning about non-RFC-compliant data.
                self.decompressor = ZLibDecompressor(
                    encoding=self.encoding, suppress_deflate_header=True
                )
            self._started_decoding = True

        low_water = self.out._low_water
        max_length = (
            0 if low_water >= sys.maxsize else max(self._max_decompress_size, low_water)
        )
        try:
            chunk = self.decompressor.decompress_sync(chunk, max_length=max_length)
        except Exception:
            raise ContentEncodingError(
                "Can not decode content-encoding: %s" % self.encoding
            )

        if chunk:
            self.out.feed_data(chunk, len(chunk))
        return self.decompressor.data_available  # type: ignore[no-any-return]

    def feed_eof(self) -> None:
        chunk = self.decompressor.flush()
        # This should never contain data as we defer the call until exhausting
        # the decompression. If .flush() is returning data, this may indicate a
        # zip bomb vulnerability as it will decompress all remaining data at once.
        assert not chunk

        if self.size > 0:
            if self.encoding == "deflate" and not self.decompressor.eof:
                raise ContentEncodingError("deflate")

        self.out.feed_eof()

    def begin_http_chunk_receiving(self) -> None:
        self.out.begin_http_chunk_receiving()

    def end_http_chunk_receiving(self) -> None:
        self.out.end_http_chunk_receiving()


HttpRequestParserPy = HttpRequestParser
HttpResponseParserPy = HttpResponseParser
RawRequestMessagePy = RawRequestMessage
RawResponseMessagePy = RawResponseMessage

try:
    if not NO_EXTENSIONS:
        from ._http_parser import (  # type: ignore[import-not-found,no-redef]
            HttpRequestParser,
            HttpResponseParser,
            RawRequestMessage,
            RawResponseMessage,
        )

        HttpRequestParserC = HttpRequestParser
        HttpResponseParserC = HttpResponseParser
        RawRequestMessageC = RawRequestMessage
        RawResponseMessageC = RawResponseMessage
except ImportError:  # pragma: no cover
    pass
