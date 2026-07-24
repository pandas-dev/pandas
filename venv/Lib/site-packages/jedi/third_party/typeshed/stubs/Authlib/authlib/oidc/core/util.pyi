from _typeshed import ReadableBuffer
from collections.abc import Iterable
from typing import SupportsBytes, SupportsIndex

def create_half_hash(
    s: str | bytes | float | Iterable[SupportsIndex] | SupportsIndex | SupportsBytes | ReadableBuffer, alg: str
) -> bytes | None: ...
