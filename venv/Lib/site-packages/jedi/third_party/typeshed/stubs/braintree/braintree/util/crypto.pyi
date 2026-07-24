from _typeshed import ReadableBuffer
from collections.abc import Iterable
from typing import Literal, overload

text_type = str

class Crypto:
    @staticmethod
    def sha1_hmac_hash(secret_key: str | ReadableBuffer, content: str | ReadableBuffer | None) -> str: ...
    @staticmethod
    def sha256_hmac_hash(secret_key: str | ReadableBuffer, content: str | ReadableBuffer | None) -> str: ...
    @overload
    @staticmethod
    def secure_compare(left: None, right: Iterable[str | bytes | bytearray]) -> Literal[False]: ...
    @overload
    @staticmethod
    def secure_compare(left: Iterable[str | bytes | bytearray], right: None) -> Literal[False]: ...
    @overload
    @staticmethod
    def secure_compare(left: None, right: None) -> Literal[False]: ...
    @overload
    @staticmethod
    def secure_compare(left: Iterable[str | bytes | bytearray], right: Iterable[str | bytes | bytearray]) -> bool: ...
