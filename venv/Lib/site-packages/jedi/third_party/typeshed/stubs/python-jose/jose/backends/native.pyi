from _typeshed import ReadableBuffer
from collections.abc import Callable
from hashlib import _Hash
from typing import Any

from .base import Key

def get_random_bytes(num_bytes: int) -> bytes: ...

class HMACKey(Key):
    HASHES: dict[str, Callable[[bytes], _Hash]]
    prepared_key: bytes
    def __init__(
        self,
        # explicitly checks for key_data as dict instance, instead of a Mapping
        key: str | bytes | dict[str, Any],
        algorithm: str,
    ) -> None: ...
    def sign(self, msg: ReadableBuffer | None) -> bytes: ...
    def verify(self, msg: ReadableBuffer | None, sig: str | bytes) -> bool: ...
    def to_dict(self) -> dict[str, Any]: ...
