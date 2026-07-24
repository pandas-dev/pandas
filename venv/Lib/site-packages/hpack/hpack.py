"""
Implements the HPACK header compression algorithm as detailed by RFC 7541.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

from .exceptions import HPACKDecodingError, InvalidTableSizeError, OversizedHeaderListError
from .huffman import HuffmanEncoder
from .huffman_constants import REQUEST_CODES, REQUEST_CODES_LENGTH
from .huffman_table import decode_huffman
from .struct import HeaderTuple, HeaderWeaklyTyped, NeverIndexedHeaderTuple
from .table import HeaderTable, table_entry_size

if TYPE_CHECKING:
    from collections.abc import Iterable  # pragma: no cover

log = logging.getLogger(__name__)

INDEX_NONE = b"\x00"
INDEX_NEVER = b"\x10"
INDEX_INCREMENTAL = b"\x40"

VARINT_MAX_LENGTH = 5 # octets, enough for encoding prefix + uint32

# Precompute 2^i for 1-8 for use in prefix calcs.
# Zero index is not used but there to save a subtraction
# as prefix numbers are not zero indexed.
_PREFIX_BIT_MAX_NUMBERS = [(2 ** i) - 1 for i in range(9)]

# We default the maximum header list we're willing to accept to 64kB. That's a
# lot of headers, but if applications want to raise it they can do.
DEFAULT_MAX_HEADER_LIST_SIZE = 2 ** 16


def _unicode_if_needed(header: HeaderWeaklyTyped, raw: bool) -> HeaderTuple:
    """
    Provides a header as a unicode string if raw is False, otherwise returns
    it as a bytestring.
    """
    name = bytes(header[0])  # type: ignore
    value = bytes(header[1])  # type: ignore

    if not raw:
        return header.__class__(name.decode("utf-8"), value.decode("utf-8"))  # type: ignore
    return header.__class__(name, value)  # type: ignore


def encode_integer(integer: int, prefix_bits: int) -> bytearray:
    """
    Encodes an integer according to the wacky integer encoding rules
    defined in the HPACK spec.
    """
    log.debug("Encoding %d with %d bits", integer, prefix_bits)

    if integer < 0:
        msg = f"Can only encode positive integers, got {integer}"
        raise ValueError(msg)

    if prefix_bits < 1 or prefix_bits > 8:
        msg = f"Prefix bits must be between 1 and 8, got {prefix_bits}"
        raise ValueError(msg)

    max_number = _PREFIX_BIT_MAX_NUMBERS[prefix_bits]

    if integer < max_number:
        return bytearray([integer])  # Seriously?
    elements = [max_number]
    integer -= max_number

    while integer >= 128:
        elements.append((integer & 127) + 128)
        integer >>= 7

    elements.append(integer)

    return bytearray(elements)


def decode_integer(data: bytes | memoryview, prefix_bits: int) -> tuple[int, int]:
    """
    Decodes an integer according to the wacky integer encoding rules
    defined in the HPACK spec. Returns a tuple of the decoded integer and the
    number of bytes that were consumed from ``data`` in order to get that
    integer.
    """
    if prefix_bits < 1 or prefix_bits > 8:
        msg = f"Prefix bits must be between 1 and 8, got {prefix_bits}"
        raise ValueError(msg)

    max_number = _PREFIX_BIT_MAX_NUMBERS[prefix_bits]
    index = 1
    shift = 0
    mask = (0xFF >> (8 - prefix_bits))

    try:
        number = data[0] & mask
        if number == max_number:
            while True:
                next_byte = data[index]
                index += 1

                if next_byte >= 128:
                    number += (next_byte - 128) << shift
                else:
                    number += next_byte << shift
                    break
                shift += 7

                if index > VARINT_MAX_LENGTH:  # sanity check to prevent infinite loops
                    # We have consumed more than enough bytes for a typical unsigned integer.
                    # The maximum size is not defined in HPACK RFC7541 Section 5.1:
                    #   > This integer representation allows for values of indefinite size.  It
                    #   > is also possible for an encoder to send a large number of zero
                    #   > values, which can waste octets and could be used to overflow integer
                    #   > values.  Integer encodings that exceed implementation limits -- in
                    #   > value or octet length -- MUST be treated as decoding errors.
                    msg = f"Variable integer representation is too long: {data!r}"
                    raise HPACKDecodingError(msg)

    except IndexError as err:
        msg = f"Unable to decode HPACK integer representation from {data!r}"
        raise HPACKDecodingError(msg) from err

    log.debug("Decoded %d, consumed %d bytes", number, index)

    return number, index


def _dict_to_iterable(header_dict: dict[bytes | str, bytes | str]) \
        -> Iterable[tuple[bytes | str, bytes | str]]:
    """
    Converts a dictionary to an iterable of key-value tuples. This is a
    HPACK-specific function because it pulls "special-headers" out first and
    then emits them.
    """
    if not isinstance(header_dict, dict):  # pragma: no cover
        msg = f"header_dict not a dict, but {type(header_dict)}"
        raise TypeError(msg)

    keys = sorted(
        header_dict.keys(),
        key=lambda k: not _to_bytes(k).startswith(b":"),
    )
    for key in keys:
        yield key, header_dict[key]


def _to_bytes(value: bytes | str | Any) -> bytes:
    """
    Convert anything to bytes through a UTF-8 encoded string
    """
    t = type(value)
    if t is bytes:
        return value  # type: ignore
    if t is not str:
        value = str(value)
    return value.encode("utf-8")  # type: ignore


class Encoder:
    """
    An HPACK encoder object. This object takes HTTP headers and emits encoded
    HTTP/2 header blocks.
    """

    def __init__(self) -> None:
        self.header_table = HeaderTable()
        self.huffman_coder = HuffmanEncoder(
            REQUEST_CODES, REQUEST_CODES_LENGTH,
        )
        self.table_size_changes: list[int] = []

    @property
    def header_table_size(self) -> int:
        """
        Controls the size of the HPACK header table.
        """
        return self.header_table.maxsize

    @header_table_size.setter
    def header_table_size(self, value: int) -> None:
        self.header_table.maxsize = value
        if self.header_table.resized:
            self.table_size_changes.append(value)

    def encode(self,
               headers: Iterable[\
                   HeaderTuple | \
                   tuple[bytes | str, bytes | str] | \
                   tuple[bytes | str, bytes | str, bool | None]] | \
                   dict[bytes | str, bytes | str],
               huffman: bool = True) -> bytes:
        """
        Takes a set of headers and encodes them into a HPACK-encoded header
        block.

        :param headers: The headers to encode. Must be either an iterable of
                        tuples, an iterable of :class:`HeaderTuple
                        <hpack.HeaderTuple>`, or a ``dict``.

                        If an iterable of tuples, the tuples may be either
                        two-tuples or three-tuples. If they are two-tuples, the
                        tuples must be of the format ``(name, value)``. If they
                        are three-tuples, they must be of the format
                        ``(name, value, sensitive)``, where ``sensitive`` is a
                        boolean value indicating whether the header should be
                        added to header tables anywhere. If not present,
                        ``sensitive`` defaults to ``False``.

                        If an iterable of :class:`HeaderTuple
                        <hpack.HeaderTuple>`, the tuples must always be
                        two-tuples. Instead of using ``sensitive`` as a third
                        tuple entry, use :class:`NeverIndexedHeaderTuple
                        <hpack.NeverIndexedHeaderTuple>` to request that
                        the field never be indexed.

                        .. warning:: HTTP/2 requires that all special headers
                            (headers whose names begin with ``:`` characters)
                            appear at the *start* of the header block. While
                            this method will ensure that happens for ``dict``
                            subclasses, callers using any other iterable of
                            tuples **must** ensure they place their special
                            headers at the start of the iterable.

                            For efficiency reasons users should prefer to use
                            iterables of two-tuples: fixing the ordering of
                            dictionary headers is an expensive operation that
                            should be avoided if possible.

        :param huffman: (optional) Whether to Huffman-encode any header sent as
                        a literal value. Except for use when debugging, it is
                        recommended that this be left enabled.

        :returns: A bytestring containing the HPACK-encoded header block.
        """
        # Transforming the headers into a header block is a procedure that can
        # be modeled as a chain or pipe. First, the headers are encoded. This
        # encoding can be done a number of ways. If the header name-value pair
        # are already in the header table we can represent them using the
        # indexed representation: the same is true if they are in the static
        # table. Otherwise, a literal representation will be used.
        header_block = []

        # Before we begin, if the header table size has been changed we need
        # to signal all changes since last emission appropriately.
        if self.header_table.resized:
            header_block.append(self._encode_table_size_change())
            self.header_table.resized = False

        if isinstance(headers, dict):
            # Turn the headers into a list of tuples if possible. This is the
            # natural way to interact with them in HPACK. Because dictionaries are
            # un-ordered, we need to make sure we grab the "special" headers first.
            hpack_headers = _dict_to_iterable(headers)
        else:
            """
            Assume headers is an iterable of HeaderTuples, or plain 2-tuples, or plain 3-tuples:

            examples:
            [
                HeaderTuple(':method', 'GET'),
                NeverIndexedHeaderTuple('customkey', 'sensitiveinfo'),
            ]
            or
            [
                (':method', 'GET'),
                ('customkey', 'some-data'),
            ]
            or
            [
                (':method', 'GET', True),
                ('customkey', 'sensitiveinfo', True),
            ]
            """
            hpack_headers = iter(headers)  # type: ignore

        # Add each header to the header block
        for header in hpack_headers:
            sensitive = False
            if isinstance(header, HeaderTuple):
                # HeaderTuple implies it's a 2-tuple with the sensitive information stored as instance attribute
                sensitive = not header.indexable
            elif len(header) > 2:
                sensitive = header[2]

            new_header = (_to_bytes(header[0]), _to_bytes(header[1]))
            header_block.append(self.add(new_header, sensitive, huffman))

        encoded = b"".join(header_block)

        log.debug("Encoded header block to %s", encoded)

        return encoded

    def add(self, to_add: tuple[bytes, bytes], sensitive: bool, huffman: bool = False) -> bytes:
        """
        Serializes a header key-value tuple.

        When sensitive is True, the header will not be added to the header table
        (see https://www.rfc-editor.org/rfc/rfc7541.html#section-7.1.3 for details),
        furthermore, the header value will be redacted in debug logs, as "SENSITIVE_REDACTED",
        to prevent accidental exposure of sensitive information.
        """
        name, value = to_add

        display_value = value if not sensitive else b"SENSITIVE_REDACTED"
        log.debug(
            "Adding %s=%s to the header table, sensitive:%s, huffman:%s",
            name,
            display_value,
            sensitive,
            huffman,
        )

        # Set our indexing mode
        indexbit = INDEX_INCREMENTAL if not sensitive else INDEX_NEVER

        # Search for a matching header in the header table.
        match = self.header_table.search(name, value)

        if match is None:
            # Not in the header table. Encode using the literal syntax,
            # and add it to the header table.
            encoded = self._encode_literal(name, value, indexbit, huffman)
            if not sensitive:
                self.header_table.add(name, value)
            return encoded

        # The header is in the table, break out the values. If we matched
        # perfectly, we can use the indexed representation: otherwise we
        # can use the indexed literal.
        index, name, perfect = match

        if perfect is not None:
            # Indexed representation.
            encoded = self._encode_indexed(index)
        else:
            # Indexed literal. We are going to add header to the
            # header table unconditionally. It is a future todo to
            # filter out headers which are known to be ineffective for
            # indexing since they just take space in the table and
            # pushed out other valuable headers.
            encoded = self._encode_indexed_literal(
                index, value, indexbit, huffman,
            )
            if not sensitive:
                self.header_table.add(name, value)

        return encoded

    def _encode_indexed(self, index: int) -> bytes:
        """
        Encodes a header using the indexed representation.
        """
        field = encode_integer(index, 7)
        field[0] |= 0x80  # we set the top bit
        return bytes(field)

    def _encode_literal(self, name: bytes, value: bytes, indexbit: bytes, huffman: bool = False) -> bytes:
        """
        Encodes a header with a literal name and literal value. If ``indexing``
        is True, the header will be added to the header table: otherwise it
        will not.
        """
        if huffman:
            name = self.huffman_coder.encode(name)
            value = self.huffman_coder.encode(value)

        name_len = encode_integer(len(name), 7)
        value_len = encode_integer(len(value), 7)

        if huffman:
            name_len[0] |= 0x80
            value_len[0] |= 0x80

        return b"".join(
            [indexbit, bytes(name_len), name, bytes(value_len), value],
        )

    def _encode_indexed_literal(self, index: int, value: bytes, indexbit: bytes, huffman: bool = False) -> bytes:
        """
        Encodes a header with an indexed name and a literal value and performs
        incremental indexing.
        """
        if indexbit != INDEX_INCREMENTAL:
            prefix = encode_integer(index, 4)
        else:
            prefix = encode_integer(index, 6)

        prefix[0] |= ord(indexbit)

        if huffman:
            value = self.huffman_coder.encode(value)

        value_len = encode_integer(len(value), 7)

        if huffman:
            value_len[0] |= 0x80

        return b"".join([bytes(prefix), bytes(value_len), value])

    def _encode_table_size_change(self) -> bytes:
        """
        Produces the encoded form of all header table size change context
        updates.
        """
        block = b""
        for size_bytes in self.table_size_changes:
            b = encode_integer(size_bytes, 5)
            b[0] |= 0x20
            block += bytes(b)
        self.table_size_changes = []
        return block


class Decoder:
    """
    An HPACK decoder object.

    .. versionchanged:: 2.3.0
       Added ``max_header_list_size`` argument.

    :param max_header_list_size: The maximum decompressed size we will allow
        for any single header block. This is a protection against DoS attacks
        that attempt to force the application to expand a relatively small
        amount of data into a really large header list, allowing enormous
        amounts of memory to be allocated.

        If this amount of data is exceeded, a `OversizedHeaderListError
        <hpack.OversizedHeaderListError>` exception will be raised. At this
        point the connection should be shut down, as the HPACK state will no
        longer be usable.

        Defaults to 64kB.
    :type max_header_list_size: ``int``
    """

    def __init__(self, max_header_list_size: int = DEFAULT_MAX_HEADER_LIST_SIZE) -> None:
        self.header_table = HeaderTable()

        #: The maximum decompressed size we will allow for any single header
        #: block. This is a protection against DoS attacks that attempt to
        #: force the application to expand a relatively small amount of data
        #: into a really large header list, allowing enormous amounts of memory
        #: to be allocated.
        #:
        #: If this amount of data is exceeded, a `OversizedHeaderListError
        #: <hpack.OversizedHeaderListError>` exception will be raised. At this
        #: point the connection should be shut down, as the HPACK state will no
        #: longer be usable.
        #:
        #: Defaults to 64kB.
        #:
        #: .. versionadded:: 2.3.0
        self.max_header_list_size = max_header_list_size

        #: Maximum allowed header table size.
        #:
        #: A HTTP/2 implementation should set this to the most recent value of
        #: SETTINGS_HEADER_TABLE_SIZE that it sent *and has received an ACK
        #: for*. Once this setting is set, the actual header table size will be
        #: checked at the end of each decoding run and whenever it is changed,
        #: to confirm that it fits in this size.
        self.max_allowed_table_size = self.header_table.maxsize

    @property
    def header_table_size(self) -> int:
        """
        Controls the size of the HPACK header table.
        """
        return self.header_table.maxsize

    @header_table_size.setter
    def header_table_size(self, value: int) -> None:
        self.header_table.maxsize = value

    def decode(self, data: bytes, raw: bool = False) -> Iterable[HeaderTuple]:
        """
        Takes an HPACK-encoded header block and decodes it into a header set.

        :param data: A bytestring representing a complete HPACK-encoded header
                     block.
        :param raw: (optional) Whether to return the headers as tuples of raw
                    byte strings or to decode them as UTF-8 before returning
                    them. The default value is False, which returns tuples of
                    Unicode strings
        :returns: A list of two-tuples of ``(name, value)`` representing the
                  HPACK-encoded headers, in the order they were decoded.
        :raises HPACKDecodingError: If an error is encountered while decoding
                                    the header block.
        """
        log.debug("Decoding %s", data)

        data_mem = memoryview(data)
        headers: list[HeaderTuple] = []
        data_len = len(data)
        inflated_size = 0
        current_index = 0

        while current_index < data_len:
            # Work out what kind of header we're decoding.
            # If the high bit is 1, it's an indexed field.
            current = data[current_index]
            indexed = bool(current & 0x80)

            # Otherwise, if the second-highest bit is 1 it's a field that does
            # alter the header table.
            literal_index = bool(current & 0x40)

            # Otherwise, if the third-highest bit is 1 it's an encoding context
            # update.
            encoding_update = bool(current & 0x20)

            if indexed:
                header, consumed = self._decode_indexed(
                    data_mem[current_index:],
                )
            elif literal_index:
                # It's a literal header that does affect the header table.
                header, consumed = self._decode_literal_index(
                    data_mem[current_index:],
                )
            elif encoding_update:
                # It's an update to the encoding context. These are forbidden
                # in a header block after any actual header.
                if headers:
                    msg = "Table size update not at the start of the block"
                    raise HPACKDecodingError(msg)
                consumed = self._update_encoding_context(
                    data_mem[current_index:],
                )
                header = None
            else:
                # It's a literal header that does not affect the header table.
                header, consumed = self._decode_literal_no_index(
                    data_mem[current_index:],
                )

            if header:
                headers.append(header)
                inflated_size += table_entry_size(header[0], header[1])

                if inflated_size > self.max_header_list_size:
                    msg = f"A header list larger than {self.max_header_list_size} has been received"
                    raise OversizedHeaderListError(msg)

            current_index += consumed

        # Confirm that the table size is lower than the maximum. We do this
        # here to ensure that we catch when the max has been *shrunk* and the
        # remote peer hasn't actually done that.
        self._assert_valid_table_size()

        try:
            return [_unicode_if_needed(h, raw) for h in headers]
        except UnicodeDecodeError as err:
            msg = "Unable to decode headers as UTF-8"
            raise HPACKDecodingError(msg) from err

    def _assert_valid_table_size(self) -> None:
        """
        Check that the table size set by the encoder is lower than the maximum
        we expect to have.
        """
        if self.header_table_size > self.max_allowed_table_size:
            msg = "Encoder did not shrink table size to within the max"
            raise InvalidTableSizeError(msg)

    def _update_encoding_context(self, data: bytes | memoryview) -> int:
        """
        Handles a byte that updates the encoding context.
        """
        # We've been asked to resize the header table.
        new_size, consumed = decode_integer(data, 5)
        if new_size > self.max_allowed_table_size:
            msg = "Encoder exceeded max allowable table size"
            raise InvalidTableSizeError(msg)
        self.header_table_size = new_size
        return consumed

    def _decode_indexed(self, data: bytes | memoryview) -> tuple[HeaderTuple, int]:
        """
        Decodes a header represented using the indexed representation.
        """
        index, consumed = decode_integer(data, 7)
        header = HeaderTuple(*self.header_table.get_by_index(index))
        log.debug("Decoded %s, consumed %d", header, consumed)
        return header, consumed

    def _decode_literal_no_index(self, data: bytes | memoryview) -> tuple[HeaderTuple, int]:
        return self._decode_literal(data, should_index=False)

    def _decode_literal_index(self, data: bytes | memoryview) -> tuple[HeaderTuple, int]:
        return self._decode_literal(data, should_index=True)

    def _decode_literal(self, data: bytes | memoryview, should_index: bool) -> tuple[HeaderTuple, int]:
        """
        Decodes a header represented with a literal.
        """
        if isinstance(data, memoryview):
            data = data.tobytes()  # pragma: no cover

        total_consumed = 0

        # When should_index is true, if the low six bits of the first byte are
        # nonzero, the header name is indexed.
        # When should_index is false, if the low four bits of the first byte
        # are nonzero the header name is indexed.
        if should_index:
            indexed_name = data[0] & 0x3F
            name_len = 6
            not_indexable = False
        else:
            high_byte = data[0]
            indexed_name = high_byte & 0x0F
            name_len = 4
            not_indexable = bool(high_byte & 0x10)

        if indexed_name:
            # Indexed header name.
            index, consumed = decode_integer(data, name_len)
            name = self.header_table.get_by_index(index)[0]

            total_consumed = consumed
            length = 0
        else:
            # Literal header name. The first byte was consumed, so we need to
            # move forward.
            data = data[1:]

            length, consumed = decode_integer(data, 7)
            name = data[consumed:consumed + length]
            if len(name) != length:
                msg = "Truncated header block"
                raise HPACKDecodingError(msg)

            if data[0] & 0x80:
                name = decode_huffman(name)
            total_consumed = consumed + length + 1  # Since we moved forward 1.

        data = data[consumed + length:]

        # The header value is definitely length-based.
        length, consumed = decode_integer(data, 7)
        value = data[consumed:consumed + length]
        if len(value) != length:
            msg = "Truncated header block"
            raise HPACKDecodingError(msg)

        if data[0] & 0x80:
            value = decode_huffman(value)

        # Updated the total consumed length.
        total_consumed += length + consumed

        # If we have been told never to index the header field, encode that in
        # the tuple we use.
        header: HeaderTuple
        if not_indexable:
            header = NeverIndexedHeaderTuple(name, value)
        else:
            header = HeaderTuple(name, value)

        # If we've been asked to index this, add it to the header table.
        if should_index:
            self.header_table.add(name, value)

        log.debug(
            "Decoded %s, total consumed %d bytes, indexed %s",
            header,
            total_consumed,
            should_index,
        )

        return header, total_consumed
