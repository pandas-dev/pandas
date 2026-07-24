from __future__ import annotations

from binascii import crc32
from collections.abc import Iterator, Mapping
from datetime import datetime
from struct import pack
from typing import TYPE_CHECKING, Any

from moto.core.model import StructureShape

if TYPE_CHECKING:
    from moto.core.model import Shape
    from moto.core.serialize import BaseEventStreamSerializer


class EncodeUtils:
    """Packing utility functions used in the encoder."""

    UINT8_BYTE_FORMAT = "!B"
    UINT16_BYTE_FORMAT = "!H"
    UINT32_BYTE_FORMAT = "!I"
    INT8_BYTE_FORMAT = "!b"
    INT16_BYTE_FORMAT = "!h"
    INT32_BYTE_FORMAT = "!i"
    INT64_BYTE_FORMAT = "!q"
    PRELUDE_BYTE_FORMAT = "!II"

    @staticmethod
    def pack_uint8(data: int) -> bytes:
        packed = pack(EncodeUtils.UINT8_BYTE_FORMAT, data)
        return packed

    @staticmethod
    def pack_uint16(data: int) -> bytes:
        packed = pack(EncodeUtils.UINT16_BYTE_FORMAT, data)
        return packed

    @staticmethod
    def pack_int32(data: int) -> bytes:
        packed = pack(EncodeUtils.INT32_BYTE_FORMAT, data)
        return packed

    @staticmethod
    def pack_int64(data: int) -> bytes:
        packed = pack(EncodeUtils.INT64_BYTE_FORMAT, data)
        return packed

    @staticmethod
    def pack_uint32(data: int) -> bytes:
        packed = pack(EncodeUtils.UINT32_BYTE_FORMAT, data)
        return packed

    @staticmethod
    def pack_byte_array(data: bytes) -> bytes:
        packed = EncodeUtils.pack_uint16(len(data)) + data
        return packed

    @staticmethod
    def pack_utf8_string(data: str) -> bytes:
        utf8_string = data.encode("utf-8")
        packed = EncodeUtils.pack_uint16(len(utf8_string)) + utf8_string
        return packed

    @staticmethod
    def pack_prelude(total_length: int, header_length: int) -> bytes:
        packed = pack(EncodeUtils.PRELUDE_BYTE_FORMAT, total_length, header_length)
        return packed


def _calculate_checksum(data: bytes, crc: int = 0) -> int:
    return crc32(data, crc) & 0xFFFFFFFF


class EventStreamEncoder:
    """Encodes Amazon event-stream binary messages.

    Wire format per
    https://smithy.io/2.0/aws/amazon-eventstream.html#amazon-event-stream-specification:

        [total_length: u32][headers_length: u32][prelude_crc: u32]
            [name_len: u8][name][type: u8][value...]    (headers, repeated)
            [payload bytes]
        [message_crc: u32]

    Both CRCs are CRC32. The message CRC is computed over the full message up
    to (but not including) the message CRC itself.
    """

    PRELUDE_LENGTH = 12
    MESSAGE_CRC_LENGTH = 4

    HEADER_TYPE_BOOL_TRUE = b"\x00"
    HEADER_TYPE_BOOL_FALSE = b"\x01"
    HEADER_TYPE_INT8 = b"\x02"
    HEADER_TYPE_INT16 = b"\x03"
    HEADER_TYPE_INT32 = b"\x04"
    HEADER_TYPE_INT64 = b"\x05"
    HEADER_TYPE_BYTE_ARRAY = b"\x06"
    HEADER_TYPE_STRING = b"\x07"
    HEADER_TYPE_TIMESTAMP = b"\x08"
    HEADER_TYPE_UUID = b"\x09"

    def encode_event(self, event: EventStreamMessage) -> bytes:
        return self.encode_message(event.headers, event.payload)

    def encode_message(self, headers: Mapping[str, Any], payload: bytes) -> bytes:
        headers_bytes = self._encode_headers(headers)
        prelude = self._encode_prelude(headers_bytes, payload)
        prelude_crc = _calculate_checksum(prelude)
        prelude_with_crc = prelude + EncodeUtils.pack_uint32(prelude_crc)
        message_without_crc = prelude_with_crc + headers_bytes + payload
        message_crc = _calculate_checksum(message_without_crc)
        return message_without_crc + EncodeUtils.pack_uint32(message_crc)

    def _encode_headers(self, headers: Mapping[str, Any]) -> bytes:
        encoded = b""
        for key, val in headers.items():
            encoded += self._encode_header_key(key)
            encoded += self._encode_header_value(val)
        return encoded

    def _encode_header_key(self, key: str) -> bytes:
        encoded_key = key.encode("utf-8")
        if len(encoded_key) > 0xFF:
            raise ValueError(f"Event-stream header name too long: {encoded_key!r}")
        return EncodeUtils.pack_uint8(len(encoded_key)) + encoded_key

    def _encode_header_value(self, value: Any) -> bytes:
        if value is True:
            return self.HEADER_TYPE_BOOL_TRUE
        elif value is False:
            return self.HEADER_TYPE_BOOL_FALSE

        if isinstance(value, int):
            return self.HEADER_TYPE_INT32 + EncodeUtils.pack_int32(value)
        elif isinstance(value, str):
            return self.HEADER_TYPE_STRING + EncodeUtils.pack_utf8_string(value)
        elif isinstance(value, (bytes, bytearray)):
            return self.HEADER_TYPE_BYTE_ARRAY + EncodeUtils.pack_byte_array(
                bytes(value)
            )
        elif isinstance(value, datetime):
            ts_ms = int(value.timestamp() * 1000)
            return self.HEADER_TYPE_TIMESTAMP + EncodeUtils.pack_int64(ts_ms)
        raise TypeError(
            f"Unsupported event-stream header value type: {type(value).__name__}"
        )

    def _encode_prelude(self, encoded_headers: bytes, payload: bytes) -> bytes:
        header_length = len(encoded_headers)
        payload_length = len(payload)
        total_length = (
            header_length
            + payload_length
            + self.PRELUDE_LENGTH
            + self.MESSAGE_CRC_LENGTH
        )
        return EncodeUtils.pack_prelude(total_length, header_length)


MESSAGE_HEADER_VALUE = bool | bytes | int | str | datetime
MESSAGE_HEADERS_DICT = dict[str, MESSAGE_HEADER_VALUE]


class EventStreamMessage:
    def __init__(self, headers: MESSAGE_HEADERS_DICT, payload: bytes) -> None:
        self.headers = headers
        self.payload = payload


class EventStream:
    """Builds an iterator of Amazon Event Stream messages from an operation result.

    The stream uses the operation output shape to locate the modeled event-stream
    member, converts each single variant event into an `EventStreamMessage`, and
    delegates payload/body serialization to the provided event-stream serializer.
    """

    def __init__(
        self,
        result: Any,
        output_shape: StructureShape,
        serializer: BaseEventStreamSerializer,
    ) -> None:
        assert output_shape.event_stream_name is not None
        self._result: Any = result
        self._output_shape: StructureShape = output_shape
        self._eventstream_member_name: str = str(output_shape.event_stream_name)
        assert isinstance(
            output_shape.members[self._eventstream_member_name], StructureShape
        )
        self._eventstream_member_shape: StructureShape = output_shape.members[
            output_shape.event_stream_name
        ]
        self._serializer: BaseEventStreamSerializer = serializer

    def __iter__(self) -> Iterator[EventStreamMessage]:
        if self._serializer.initial_response_required:
            yield self.get_initial_response()
        event_stream = self._serializer.get_value(
            self._result, self._eventstream_member_name, self._eventstream_member_shape
        )
        if event_stream is not None:
            for event in event_stream:
                if not isinstance(event, Mapping) or len(event) != 1:
                    raise ValueError(
                        "Each event-stream item must be a single-key mapping "
                        "of {variant_name: value}"
                    )
                variant_name, variant_value = next(iter(event.items()))
                variant_shape = self._eventstream_member_shape.members[variant_name]
                assert isinstance(variant_shape, StructureShape)
                yield self._generate_event(variant_name, variant_value, variant_shape)

    def _generate_event(
        self, variant_name: str, value: Any, variant_shape: StructureShape
    ) -> EventStreamMessage:
        headers: MESSAGE_HEADERS_DICT = {
            ":message-type": "event",
            ":event-type": variant_name,
        }
        eventpayload_member: tuple[str, Shape] | None = None
        body_members: list[tuple[str, Shape]] = []
        for mn, ms in variant_shape.members.items():
            if ms.serialization.get("eventpayload"):
                eventpayload_member = (mn, ms)
            else:
                body_members.append((mn, ms))
        payload: bytes = b""
        content_type: str | None = None
        if eventpayload_member is not None:
            mn, ms = eventpayload_member
            member_value = self._serializer.get_value(value, mn, ms)
            if ms.type_name == "blob":
                if member_value is None:
                    payload = b""
                elif isinstance(member_value, (bytes, bytearray)):
                    payload = bytes(member_value)
                else:
                    payload = str(member_value).encode("utf-8")
                content_type = "application/octet-stream"
            elif ms.type_name == "string":
                payload = (member_value or "").encode("utf-8")
                content_type = "text/plain"
            else:
                assert isinstance(ms, StructureShape)
                payload, content_type = self._serializer.serialize_event_payload(
                    member_value, ms, list(ms.members.items())
                )
        elif body_members:
            payload, content_type = self._serializer.serialize_event_payload(
                value, variant_shape, body_members
            )
        if content_type:
            headers[":content-type"] = content_type
        return EventStreamMessage(headers, payload)

    def get_initial_response(self) -> EventStreamMessage:
        # Initial-response covers any output members not part of the event stream.
        initial_response_members = [
            (mn, ms)
            for mn, ms in self._output_shape.members.items()
            if mn != self._eventstream_member_name
        ]
        payload, content_type = b"", None
        if initial_response_members:
            payload, content_type = self._serializer.serialize_event_payload(
                self._result, self._output_shape, initial_response_members
            )
        headers: MESSAGE_HEADERS_DICT = {
            ":message-type": "event",
            ":event-type": "initial-response",
        }
        if content_type:
            headers[":content-type"] = content_type
        return EventStreamMessage(headers, payload)
