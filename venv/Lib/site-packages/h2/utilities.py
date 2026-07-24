"""
h2/utilities
~~~~~~~~~~~~

Utility functions that do not belong in a separate module.
"""
from __future__ import annotations

import collections
from typing import TYPE_CHECKING, Any, NamedTuple

from hpack.struct import HeaderTuple, NeverIndexedHeaderTuple

from .exceptions import FlowControlError, ProtocolError

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Generator, Iterable

    from hpack.struct import Header, HeaderWeaklyTyped

SIGIL = ord(b":")
INFORMATIONAL_START = ord(b"1")


# A set of headers that are hop-by-hop or connection-specific and thus
# forbidden in HTTP/2. This list comes from RFC 7540 § 8.1.2.2.
CONNECTION_HEADERS = frozenset([
    b"connection",
    b"proxy-connection",
    b"keep-alive",
    b"transfer-encoding",
    b"upgrade",
])


_ALLOWED_PSEUDO_HEADER_FIELDS = frozenset([
    b":method",
    b":scheme",
    b":authority",
    b":path",
    b":status",
    b":protocol",
])


_SECURE_HEADERS = frozenset([
    # May have basic credentials which are vulnerable to dictionary attacks.
    b"authorization",
    b"proxy-authorization",
])


_REQUEST_ONLY_HEADERS = frozenset([
    b":scheme",
    b":path",
    b":authority",
    b":method",
    b":protocol",
])


_RESPONSE_ONLY_HEADERS = frozenset([b":status"])


# A Set of pseudo headers that are only valid if the method is
# CONNECT, see RFC 8441 § 5
_CONNECT_REQUEST_ONLY_HEADERS = frozenset([b":protocol"])


def _secure_headers(headers: Iterable[Header],
                    hdr_validation_flags: HeaderValidationFlags | None) -> Generator[Header, None, None]:
    """
    Certain headers are at risk of being attacked during the header compression
    phase, and so need to be kept out of header compression contexts. This
    function automatically transforms certain specific headers into HPACK
    never-indexed fields to ensure they don't get added to header compression
    contexts.

    This function currently implements two rules:

    - 'authorization' and 'proxy-authorization' fields are automatically made
      never-indexed.
    - Any 'cookie' header field shorter than 20 bytes long is made
      never-indexed.

    These fields are the most at-risk. These rules are inspired by Firefox
    and nghttp2.
    """
    for header in headers:
        assert isinstance(header[0], bytes)
        if header[0] in _SECURE_HEADERS or (header[0] in b"cookie" and len(header[1]) < 20):
            yield NeverIndexedHeaderTuple(header[0], header[1])
        else:
            yield header


def extract_method_header(headers: Iterable[Header]) -> bytes | None:
    """
    Extracts the request method from the headers list.
    """
    for k, v in headers:
        if isinstance(v, bytes) and k == b":method":
            return v
        if isinstance(v, str) and k == ":method":
            return v.encode("utf-8")  # pragma: no cover
    return None


def is_informational_response(headers: Iterable[Header]) -> bool:
    """
    Searches headers list for a :status header to confirm that a given
    collection of headers are an informational response. Assumes the header
    are well formed and encoded as bytes: that is, that the HTTP/2 special
    headers are first in the block, and so that it can stop looking when it
    finds the first header field whose name does not begin with a colon.

    :param headers: The HTTP/2 headers.
    :returns: A boolean indicating if this is an informational response.
    """
    for n, v in headers:
        if not n.startswith(b":"):
            return False
        if n != b":status":
            # If we find a non-special header, we're done here: stop looping.
            continue
        # If the first digit is a 1, we've got informational headers.
        return v.startswith(b"1")
    return False


def guard_increment_window(current: int, increment: int) -> int:
    """
    Increments a flow control window, guarding against that window becoming too
    large.

    :param current: The current value of the flow control window.
    :param increment: The increment to apply to that window.
    :returns: The new value of the window.
    :raises: ``FlowControlError``
    """
    # The largest value the flow control window may take.
    LARGEST_FLOW_CONTROL_WINDOW = 2**31 - 1  # noqa: N806

    new_size = current + increment

    if new_size > LARGEST_FLOW_CONTROL_WINDOW:
        msg = f"May not increment flow control window past {LARGEST_FLOW_CONTROL_WINDOW}"
        raise FlowControlError(msg)

    return new_size


def authority_from_headers(headers: Iterable[Header]) -> bytes | None:
    """
    Given a header set, searches for the authority header and returns the
    value.

    Note that this doesn't use indexing, so should only be called if the
    headers are for a client request. Otherwise, will loop over the entire
    header set, which is potentially unwise.

    :param headers: The HTTP header set.
    :returns: The value of the authority header, or ``None``.
    :rtype: ``bytes`` or ``None``.
    """
    for n, v in headers:
        if n == b":authority":
            return v

    return None


# Flags used by the validate_headers pipeline to determine which checks
# should be applied to a given set of headers.
class HeaderValidationFlags(NamedTuple):
    is_client: bool
    is_trailer: bool
    is_response_header: bool
    is_push_promise: bool


def validate_headers(headers: Iterable[Header], hdr_validation_flags: HeaderValidationFlags) -> Iterable[Header]:
    """
    Validates a header sequence against a set of constraints from RFC 7540.

    :param headers: The HTTP header set.
    :param hdr_validation_flags: An instance of HeaderValidationFlags.
    """
    # This validation logic is built on a sequence of generators that are
    # iterated over to provide the final header list. This reduces some of the
    # overhead of doing this checking. However, it's worth noting that this
    # checking remains somewhat expensive, and attempts should be made wherever
    # possible to reduce the time spent doing them.
    #
    # For example, we avoid tuple unpacking in loops because it represents a
    # fixed cost that we don't want to spend, instead indexing into the header
    # tuples.
    headers = _reject_illegal_characters(
        headers, hdr_validation_flags,
    )
    headers = _reject_empty_header_names(
        headers, hdr_validation_flags,
    )
    headers = _reject_te(
        headers, hdr_validation_flags,
    )
    headers = _reject_connection_header(
        headers, hdr_validation_flags,
    )
    headers = _reject_pseudo_header_fields(
        headers, hdr_validation_flags,
    )
    headers = _check_host_authority_header(
        headers, hdr_validation_flags,
    )
    return _check_path_header(headers, hdr_validation_flags)


def _reject_illegal_characters(headers: Iterable[Header],
                               hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    Raises a ProtocolError if any header names or values contain illegal characters.
    See <https://www.rfc-editor.org/rfc/rfc9113.html#section-8.2.1>.
    """
    for header in headers:
        # > A field name MUST NOT contain characters in the ranges 0x00-0x20, 0x41-0x5a,
        # > or 0x7f-0xff (all ranges inclusive).
        for c in header[0]:
            if 0x41 <= c <= 0x5a:
                msg = f"Received uppercase header name {header[0]!r}."
                raise ProtocolError(msg)
            if c <= 0x20 or c >= 0x7f:
                msg = f"Illegal character '{chr(c)}' in header name: {header[0]!r}"
                raise ProtocolError(msg)

        # > With the exception of pseudo-header fields (Section 8.3), which have a name
        # > that starts with a single colon, field names MUST NOT include a colon (ASCII
        # > COLON, 0x3a).
        if header[0].find(b":", 1) != -1:
            msg = f"Illegal character ':' in header name: {header[0]!r}"
            raise ProtocolError(msg)

        # For compatibility with RFC 7230 header fields, we need to allow the field
        # value to be an empty string. This is ludicrous, but technically allowed.
        if field_value := header[1]:

            # > A field value MUST NOT contain the zero value (ASCII NUL, 0x00), line feed
            # > (ASCII LF, 0x0a), or carriage return (ASCII CR, 0x0d) at any position.
            for c in field_value:
                if c == 0 or c == 0x0a or c == 0x0d:  # noqa: PLR1714
                    msg = f"Illegal character '{chr(c)}' in header value: {field_value!r}"
                    raise ProtocolError(msg)

            # > A field value MUST NOT start or end with an ASCII whitespace character
            # > (ASCII SP or HTAB, 0x20 or 0x09).
            if (
                field_value[0] == 0x20 or
                field_value[0] == 0x09 or
                field_value[-1] == 0x20 or
                field_value[-1] == 0x09
            ):
                msg = f"Received header value surrounded by whitespace {field_value!r}"
                raise ProtocolError(msg)

        yield header


def _reject_empty_header_names(headers: Iterable[Header],
                               hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    Raises a ProtocolError if any header names are empty (length 0).
    While hpack decodes such headers without errors, they are semantically
    forbidden in HTTP, see RFC 7230, stating that they must be at least one
    character long.
    """
    for header in headers:
        if len(header[0]) == 0:
            msg = "Received header name with zero length."
            raise ProtocolError(msg)
        yield header


def _reject_te(headers: Iterable[Header], hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    Raises a ProtocolError if the TE header is present in a header block and
    its value is anything other than "trailers".
    """
    for header in headers:
        if header[0] == b"te" and header[1].lower() != b"trailers":
            msg = f"Invalid value for TE header: {header[1]!r}"
            raise ProtocolError(msg)

        yield header


def _reject_connection_header(headers: Iterable[Header], hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    Raises a ProtocolError if the Connection header is present in a header
    block.
    """
    for header in headers:
        if header[0] in CONNECTION_HEADERS:
            msg = f"Connection-specific header field present: {header[0]!r}."
            raise ProtocolError(msg)

        yield header


def _assert_header_in_set(bytes_header: bytes,
                          header_set: set[bytes | str] | set[bytes] | set[str]) -> None:
    """
    Given a set of header names, checks whether the string or byte version of
    the header name is present. Raises a Protocol error with the appropriate
    error if it's missing.
    """
    if bytes_header not in header_set:
        msg = f"Header block missing mandatory {bytes_header!r} header"
        raise ProtocolError(msg)


def _reject_pseudo_header_fields(headers: Iterable[Header],
                                 hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    Raises a ProtocolError if duplicate pseudo-header fields are found in a
    header block or if a pseudo-header field appears in a block after an
    ordinary header field.

    Raises a ProtocolError if pseudo-header fields are found in trailers.
    """
    seen_pseudo_header_fields = set()
    seen_regular_header = False
    method = None

    for header in headers:
        if header[0][0] == SIGIL:
            if header[0] in seen_pseudo_header_fields:
                msg = f"Received duplicate pseudo-header field {header[0]!r}"
                raise ProtocolError(msg)

            seen_pseudo_header_fields.add(header[0])

            if seen_regular_header:
                msg = f"Received pseudo-header field out of sequence: {header[0]!r}"
                raise ProtocolError(msg)

            if header[0] not in _ALLOWED_PSEUDO_HEADER_FIELDS:
                msg = f"Received custom pseudo-header field {header[0]!r}"
                raise ProtocolError(msg)

            if header[0] in b":method":
                method = header[1]

        else:
            seen_regular_header = True

        yield header

    # Check the pseudo-headers we got to confirm they're acceptable.
    _check_pseudo_header_field_acceptability(
        seen_pseudo_header_fields, method, hdr_validation_flags,
    )


def _check_pseudo_header_field_acceptability(pseudo_headers: set[bytes | str] | set[bytes] | set[str],
                                             method: bytes | None,
                                             hdr_validation_flags: HeaderValidationFlags) -> None:
    """
    Given the set of pseudo-headers present in a header block and the
    validation flags, confirms that RFC 7540 allows them.
    """
    # Pseudo-header fields MUST NOT appear in trailers - RFC 7540 § 8.1.2.1
    if hdr_validation_flags.is_trailer and pseudo_headers:
        msg = f"Received pseudo-header in trailer {pseudo_headers}"
        raise ProtocolError(msg)

    # If ':status' pseudo-header is not there in a response header, reject it.
    # Similarly, if ':path', ':method', or ':scheme' are not there in a request
    # header, reject it. Additionally, if a response contains any request-only
    # headers or vice-versa, reject it.
    # Relevant RFC section: RFC 7540 § 8.1.2.4
    # https://tools.ietf.org/html/rfc7540#section-8.1.2.4
    if hdr_validation_flags.is_response_header:
        _assert_header_in_set(b":status", pseudo_headers)
        invalid_response_headers = pseudo_headers & _REQUEST_ONLY_HEADERS
        if invalid_response_headers:
            msg = f"Encountered request-only headers {invalid_response_headers}"
            raise ProtocolError(msg)
    elif (not hdr_validation_flags.is_response_header and
          not hdr_validation_flags.is_trailer):
        # This is a request, so we need to have seen :path, :method, and
        # :scheme.
        _assert_header_in_set(b":path", pseudo_headers)
        _assert_header_in_set(b":method", pseudo_headers)
        _assert_header_in_set(b":scheme", pseudo_headers)
        invalid_request_headers = pseudo_headers & _RESPONSE_ONLY_HEADERS
        if invalid_request_headers:
            msg = f"Encountered response-only headers {invalid_request_headers}"
            raise ProtocolError(msg)
        if method != b"CONNECT":
            invalid_headers = pseudo_headers & _CONNECT_REQUEST_ONLY_HEADERS
            if invalid_headers:
                msg = f"Encountered connect-request-only headers {invalid_headers!r}"
                raise ProtocolError(msg)


def _validate_host_authority_header(headers: Iterable[Header]) -> Generator[Header, None, None]:
    """
    Given the :authority and Host headers from a request block that isn't
    a trailer, check that:
     1. At least one of these headers is set.
     2. If both headers are set, they match.

    :param headers: The HTTP header set.
    :raises: ``ProtocolError``
    """
    # We use None as a sentinel value.  Iterate over the list of headers,
    # and record the value of these headers (if present).  We don't need
    # to worry about receiving duplicate :authority headers, as this is
    # enforced by the _reject_pseudo_header_fields() pipeline.
    #
    # TODO: We should also guard against receiving duplicate Host headers,
    # and against sending duplicate headers.
    authority_header_val = None
    host_header_val = None

    for header in headers:
        if header[0] == b":authority":
            authority_header_val = header[1]
        elif header[0] == b"host":
            host_header_val = header[1]

        yield header

    # If we have not-None values for these variables, then we know we saw
    # the corresponding header.
    authority_present = (authority_header_val is not None)
    host_present = (host_header_val is not None)

    # It is an error for a request header block to contain neither
    # an :authority header nor a Host header.
    if not authority_present and not host_present:
        msg = "Request header block does not have an :authority or Host header."
        raise ProtocolError(msg)

    # If we receive both headers, they should definitely match.
    if authority_present and host_present and authority_header_val != host_header_val:
        msg = (
            "Request header block has mismatched :authority and "
            f"Host headers: {authority_header_val!r} / {host_header_val!r}"
        )
        raise ProtocolError(msg)


def _check_host_authority_header(headers: Iterable[Header],
                                 hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    Raises a ProtocolError if a header block arrives that does not contain an
    :authority or a Host header, or if a header block contains both fields,
    but their values do not match.
    """
    # We only expect to see :authority and Host headers on request header
    # blocks that aren't trailers, so skip this validation if this is a
    # response header or we're looking at trailer blocks.
    skip_validation = (
        hdr_validation_flags.is_response_header or
        hdr_validation_flags.is_trailer
    )
    if skip_validation:
        return (h for h in headers)

    return _validate_host_authority_header(headers)


def _check_path_header(headers: Iterable[Header],
                       hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    Raise a ProtocolError if a header block arrives or is sent that contains an
    empty :path header.
    """
    def inner() -> Generator[Header, None, None]:
        for header in headers:
            if header[0] == b":path" and not header[1]:
                msg = "An empty :path header is forbidden"
                raise ProtocolError(msg)

            yield header

    # We only expect to see :authority and Host headers on request header
    # blocks that aren't trailers, so skip this validation if this is a
    # response header or we're looking at trailer blocks.
    skip_validation = (
        hdr_validation_flags.is_response_header or
        hdr_validation_flags.is_trailer
    )
    if skip_validation:
        return (h for h in headers)
    return inner()


def _to_bytes(v: bytes | str) -> bytes:
    """
    Given an assumed `str` (or anything that supports `.encode()`),
    encodes it using utf-8 into bytes. Returns the unmodified object
    if it is already a `bytes` object.
    """
    return v if isinstance(v, bytes) else v.encode("utf-8")


def utf8_encode_headers(headers: Iterable[HeaderWeaklyTyped]) -> list[Header]:
    """
    Given an iterable of header two-tuples, rebuilds that as a list with the
    header names and values encoded as utf-8 bytes. This function produces
    tuples that preserve the original type of the header tuple for tuple and
    any ``HeaderTuple``.
    """
    encoded_headers: list[Header] = []
    for header in headers:
        h = (_to_bytes(header[0]), _to_bytes(header[1]))
        if isinstance(header, HeaderTuple):
            encoded_headers.append(header.__class__(h[0], h[1]))
        else:
            encoded_headers.append(h)
    return encoded_headers


def _lowercase_header_names(headers: Iterable[Header],
                            hdr_validation_flags: HeaderValidationFlags | None) -> Generator[Header, None, None]:
    """
    Given an iterable of header two-tuples, rebuilds that iterable with the
    header names lowercased. This generator produces tuples that preserve the
    original type of the header tuple for tuple and any ``HeaderTuple``.
    """
    for header in headers:
        if isinstance(header, HeaderTuple):
            yield header.__class__(header[0].lower(), header[1])
        else:
            yield (header[0].lower(), header[1])


def _strip_surrounding_whitespace(headers: Iterable[Header],
                                  hdr_validation_flags: HeaderValidationFlags | None) -> Generator[Header, None, None]:
    """
    Given an iterable of header two-tuples, strip both leading and trailing
    whitespace from both header names and header values. This generator
    produces tuples that preserve the original type of the header tuple for
    tuple and any ``HeaderTuple``.
    """
    for header in headers:
        if isinstance(header, HeaderTuple):
            yield header.__class__(header[0].strip(), header[1].strip())
        else:
            yield (header[0].strip(), header[1].strip())


def _strip_connection_headers(headers: Iterable[Header],
                              hdr_validation_flags: HeaderValidationFlags | None) -> Generator[Header, None, None]:
    """
    Strip any connection headers as per RFC7540 § 8.1.2.2.
    """
    for header in headers:
        if header[0] not in CONNECTION_HEADERS:
            yield header


def _check_sent_host_authority_header(headers: Iterable[Header],
                                      hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    Raises an InvalidHeaderBlockError if we try to send a header block
    that does not contain an :authority or a Host header, or if
    the header block contains both fields, but their values do not match.
    """
    # We only expect to see :authority and Host headers on request header
    # blocks that aren't trailers, so skip this validation if this is a
    # response header or we're looking at trailer blocks.
    skip_validation = (
        hdr_validation_flags.is_response_header or
        hdr_validation_flags.is_trailer
    )
    if skip_validation:
        return (h for h in headers)

    return _validate_host_authority_header(headers)


def _combine_cookie_fields(headers: Iterable[Header], hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    RFC 7540 § 8.1.2.5 allows HTTP/2 clients to split the Cookie header field,
    which must normally appear only once, into multiple fields for better
    compression. However, they MUST be joined back up again when received.
    This normalization step applies that transform. The side-effect is that
    all cookie fields now appear *last* in the header block.
    """
    # There is a problem here about header indexing. Specifically, it's
    # possible that all these cookies are sent with different header indexing
    # values. At this point it shouldn't matter too much, so we apply our own
    # logic and make them never-indexed.
    cookies: list[bytes] = []
    for header in headers:
        if header[0] == b"cookie":
            cookies.append(header[1])
        else:
            yield header
    if cookies:
        cookie_val = b"; ".join(cookies)
        yield NeverIndexedHeaderTuple(b"cookie", cookie_val)


def _split_outbound_cookie_fields(headers: Iterable[Header],
                                  hdr_validation_flags: HeaderValidationFlags | None) -> Generator[Header, None, None]:
    """
    RFC 7540 § 8.1.2.5 allows for better compression efficiency,
    to split the Cookie header field into separate header fields

    We want to do it for outbound requests, as we are doing for
    inbound.
    """
    for header in headers:
        assert isinstance(header[0], bytes)
        assert isinstance(header[1], bytes)
        if header[0] == b"cookie":
            for cookie_val in header[1].split(b"; "):
                if isinstance(header, HeaderTuple):
                    yield header.__class__(header[0], cookie_val)
                else:
                    yield header[0], cookie_val
        else:
            yield header


def normalize_outbound_headers(headers: Iterable[Header],
                               hdr_validation_flags: HeaderValidationFlags | None,
                               should_split_outbound_cookies: bool=False) -> Generator[Header, None, None]:
    """
    Normalizes a header sequence that we are about to send.

    :param headers: The HTTP header set.
    :param hdr_validation_flags: An instance of HeaderValidationFlags.
    :param should_split_outbound_cookies: boolean flag
    """
    headers = _lowercase_header_names(headers, hdr_validation_flags)
    if should_split_outbound_cookies:
        headers = _split_outbound_cookie_fields(headers, hdr_validation_flags)
    headers = _strip_surrounding_whitespace(headers, hdr_validation_flags)
    headers = _strip_connection_headers(headers, hdr_validation_flags)
    return _secure_headers(headers, hdr_validation_flags)



def normalize_inbound_headers(headers: Iterable[Header],
                              hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    Normalizes a header sequence that we have received.

    :param headers: The HTTP header set.
    :param hdr_validation_flags: An instance of HeaderValidationFlags
    """
    return _combine_cookie_fields(headers, hdr_validation_flags)


def validate_outbound_headers(headers: Iterable[Header],
                              hdr_validation_flags: HeaderValidationFlags) -> Generator[Header, None, None]:
    """
    Validates and normalizes a header sequence that we are about to send.

    :param headers: The HTTP header set.
    :param hdr_validation_flags: An instance of HeaderValidationFlags.
    """
    headers = _reject_te(
        headers, hdr_validation_flags,
    )
    headers = _reject_connection_header(
        headers, hdr_validation_flags,
    )
    headers = _reject_pseudo_header_fields(
        headers, hdr_validation_flags,
    )
    headers = _check_sent_host_authority_header(
        headers, hdr_validation_flags,
    )
    return _check_path_header(headers, hdr_validation_flags)



class SizeLimitDict(collections.OrderedDict[int, Any]):

    def __init__(self, *args: dict[int, int], **kwargs: Any) -> None:
        self._size_limit = kwargs.pop("size_limit", None)
        super().__init__(*args, **kwargs)

        self._check_size_limit()

    def __setitem__(self, key: int, value: Any | int) -> None:
        super().__setitem__(key, value)

        self._check_size_limit()

    def _check_size_limit(self) -> None:
        if self._size_limit is not None:
            while len(self) > self._size_limit:
                self.popitem(last=False)
