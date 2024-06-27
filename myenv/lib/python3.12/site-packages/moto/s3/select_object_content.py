import binascii
import struct
from typing import Any, Dict, List, Optional


def parse_query(text_input: str, query: str) -> List[Dict[str, Any]]:
    from py_partiql_parser import S3SelectParser

    return S3SelectParser(source_data={"s3object": text_input}).parse(query)


def _create_header(key: bytes, value: bytes) -> bytes:
    return struct.pack("b", len(key)) + key + struct.pack("!bh", 7, len(value)) + value


def _create_message(
    content_type: Optional[bytes], event_type: bytes, payload: bytes
) -> bytes:
    headers = _create_header(b":message-type", b"event")
    headers += _create_header(b":event-type", event_type)
    if content_type is not None:
        headers += _create_header(b":content-type", content_type)

    headers_length = struct.pack("!I", len(headers))
    total_length = struct.pack("!I", len(payload) + len(headers) + 16)
    prelude = total_length + headers_length

    prelude_crc = struct.pack("!I", binascii.crc32(total_length + headers_length))
    message_crc = struct.pack(
        "!I", binascii.crc32(prelude + prelude_crc + headers + payload)
    )

    return prelude + prelude_crc + headers + payload + message_crc


def _create_stats_message() -> bytes:
    stats = b"""<Stats><BytesScanned>24</BytesScanned><BytesProcessed>24</BytesProcessed><BytesReturned>22</BytesReturned></Stats>"""
    return _create_message(content_type=b"text/xml", event_type=b"Stats", payload=stats)


def _create_data_message(payload: bytes) -> bytes:
    # https://docs.aws.amazon.com/AmazonS3/latest/API/RESTSelectObjectAppendix.html
    return _create_message(
        content_type=b"application/octet-stream", event_type=b"Records", payload=payload
    )


def _create_end_message() -> bytes:
    return _create_message(content_type=None, event_type=b"End", payload=b"")


def serialize_select(data_list: List[bytes]) -> bytes:
    response = b""
    for data in data_list:
        response += _create_data_message(data)
    return response + _create_stats_message() + _create_end_message()
