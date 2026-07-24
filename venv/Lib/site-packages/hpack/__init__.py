"""
HTTP/2 header encoding for Python.
"""
from __future__ import annotations

from .exceptions import HPACKDecodingError, HPACKError, InvalidTableIndex, InvalidTableIndexError, InvalidTableSizeError, OversizedHeaderListError
from .hpack import Decoder, Encoder
from .struct import HeaderTuple, NeverIndexedHeaderTuple

__all__ = [
    "Decoder",
    "Encoder",
    "HPACKDecodingError",
    "HPACKError",
    "HeaderTuple",
    "InvalidTableIndex",
    "InvalidTableIndexError",
    "InvalidTableSizeError",
    "NeverIndexedHeaderTuple",
    "OversizedHeaderListError",
]

__version__ = "4.2.0"
