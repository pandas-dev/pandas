from _typeshed import Incomplete
from typing import Final

TYPE_NONE: Final = 0x00
TYPE_BOOL: Final = 0x01
TYPE_INT64: Final = 0x05
TYPE_FLOAT64: Final = 0x06
TYPE_BYTES: Final = 0x10
TYPE_STRING: Final = 0x11
TYPE_LIST: Final = 0x20
TYPE_DICT: Final = 0x21
MAX_STRING_SIZE: Final = 67108864
MAX_BYTES_SIZE: Final = 67108864
MAX_LIST_SIZE: Final = 1048576
MAX_DICT_SIZE: Final = 1048576

class TLVEncoder:
    @staticmethod
    def encode(
        value: bool | float | str | bytes | list[Incomplete] | tuple[Incomplete, ...] | dict[object, Incomplete] | None,
    ) -> bytes: ...  # dict key passed to `str()` function
    @staticmethod
    def decode(data: bytes, offset: int = 0) -> tuple[Incomplete, int]: ...
    @staticmethod
    def decode_full(data: bytes): ...
