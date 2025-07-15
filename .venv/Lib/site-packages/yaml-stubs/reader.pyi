from _typeshed import SupportsRead
from typing import Any
from typing_extensions import TypeAlias

from yaml.error import YAMLError

_ReadStream: TypeAlias = str | bytes | SupportsRead[str] | SupportsRead[bytes]

class ReaderError(YAMLError):
    name: Any
    character: Any
    position: Any
    encoding: Any
    reason: Any
    def __init__(self, name, position, character, encoding, reason) -> None: ...

class Reader:
    name: Any
    stream: SupportsRead[str] | SupportsRead[bytes] | None
    stream_pointer: Any
    eof: Any
    buffer: Any
    pointer: Any
    raw_buffer: Any
    raw_decode: Any
    encoding: Any
    index: Any
    line: Any
    column: Any
    def __init__(self, stream: _ReadStream) -> None: ...
    def peek(self, index=0): ...
    def prefix(self, length=1): ...
    def forward(self, length=1): ...
    def get_mark(self): ...
    def determine_encoding(self): ...
    NON_PRINTABLE: Any
    def check_printable(self, data): ...
    def update(self, length): ...
    def update_raw(self, size=4096): ...

__all__ = ["Reader", "ReaderError"]
