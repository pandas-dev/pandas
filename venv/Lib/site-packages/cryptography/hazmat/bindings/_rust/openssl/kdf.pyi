# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

import typing

from cryptography.hazmat.primitives.hashes import HashAlgorithm
from cryptography.hazmat.primitives.kdf.kbkdf import CounterLocation, Mode
from cryptography.utils import Buffer

class PBKDF2HMAC:
    def __init__(
        self,
        algorithm: HashAlgorithm,
        length: int,
        salt: bytes,
        iterations: int,
        backend: typing.Any = None,
    ) -> None: ...
    def derive(self, key_material: Buffer) -> bytes: ...
    def derive_into(self, key_material: Buffer, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...

class Scrypt:
    def __init__(
        self,
        salt: bytes,
        length: int,
        n: int,
        r: int,
        p: int,
        backend: typing.Any = None,
    ) -> None: ...
    def derive(self, key_material: Buffer) -> bytes: ...
    def derive_into(self, key_material: Buffer, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...

class Argon2d:
    def __init__(
        self,
        *,
        salt: bytes,
        length: int,
        iterations: int,
        lanes: int,
        memory_cost: int,
        ad: bytes | None = None,
        secret: bytes | None = None,
    ) -> None: ...
    def derive(self, key_material: bytes) -> bytes: ...
    def derive_into(self, key_material: bytes, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...
    def derive_phc_encoded(self, key_material: bytes) -> str: ...
    @classmethod
    def verify_phc_encoded(
        cls, key_material: bytes, phc_encoded: str, secret: bytes | None = None
    ) -> None: ...

class Argon2i:
    def __init__(
        self,
        *,
        salt: bytes,
        length: int,
        iterations: int,
        lanes: int,
        memory_cost: int,
        ad: bytes | None = None,
        secret: bytes | None = None,
    ) -> None: ...
    def derive(self, key_material: bytes) -> bytes: ...
    def derive_into(self, key_material: bytes, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...
    def derive_phc_encoded(self, key_material: bytes) -> str: ...
    @classmethod
    def verify_phc_encoded(
        cls, key_material: bytes, phc_encoded: str, secret: bytes | None = None
    ) -> None: ...

class Argon2id:
    def __init__(
        self,
        *,
        salt: bytes,
        length: int,
        iterations: int,
        lanes: int,
        memory_cost: int,
        ad: bytes | None = None,
        secret: bytes | None = None,
    ) -> None: ...
    def derive(self, key_material: bytes) -> bytes: ...
    def derive_into(self, key_material: bytes, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...
    def derive_phc_encoded(self, key_material: bytes) -> str: ...
    @classmethod
    def verify_phc_encoded(
        cls, key_material: bytes, phc_encoded: str, secret: bytes | None = None
    ) -> None: ...

class HKDF:
    def __init__(
        self,
        algorithm: HashAlgorithm,
        length: int,
        salt: bytes | None,
        info: bytes | None,
        backend: typing.Any = None,
    ): ...
    @staticmethod
    def extract(
        algorithm: HashAlgorithm, salt: bytes | None, key_material: Buffer
    ) -> bytes: ...
    def derive(self, key_material: Buffer) -> bytes: ...
    def derive_into(self, key_material: Buffer, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...

class HKDFExpand:
    def __init__(
        self,
        algorithm: HashAlgorithm,
        length: int,
        info: bytes | None,
        backend: typing.Any = None,
    ): ...
    def derive(self, key_material: Buffer) -> bytes: ...
    def derive_into(self, key_material: Buffer, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...

class X963KDF:
    def __init__(
        self,
        algorithm: HashAlgorithm,
        length: int,
        sharedinfo: bytes | None,
        backend: typing.Any = None,
    ) -> None: ...
    def derive(self, key_material: Buffer) -> bytes: ...
    def derive_into(self, key_material: Buffer, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...

class ConcatKDFHash:
    def __init__(
        self,
        algorithm: HashAlgorithm,
        length: int,
        otherinfo: bytes | None,
        backend: typing.Any = None,
    ) -> None: ...
    def derive(self, key_material: Buffer) -> bytes: ...
    def derive_into(self, key_material: Buffer, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...

class ConcatKDFHMAC:
    def __init__(
        self,
        algorithm: HashAlgorithm,
        length: int,
        salt: bytes | None,
        otherinfo: bytes | None,
        backend: typing.Any = None,
    ) -> None: ...
    def derive(self, key_material: Buffer) -> bytes: ...
    def derive_into(self, key_material: Buffer, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...

class KBKDFHMAC:
    def __init__(
        self,
        algorithm: HashAlgorithm,
        mode: Mode,
        length: int,
        rlen: int,
        llen: int | None,
        location: CounterLocation,
        label: bytes | None,
        context: bytes | None,
        fixed: bytes | None,
        backend: typing.Any = None,
        *,
        break_location: int | None = None,
    ) -> None: ...
    def derive(self, key_material: Buffer) -> bytes: ...
    def derive_into(self, key_material: Buffer, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...

class KBKDFCMAC:
    def __init__(
        self,
        algorithm: typing.Any,
        mode: Mode,
        length: int,
        rlen: int,
        llen: int | None,
        location: CounterLocation,
        label: bytes | None,
        context: bytes | None,
        fixed: bytes | None,
        backend: typing.Any = None,
        *,
        break_location: int | None = None,
    ) -> None: ...
    def derive(self, key_material: Buffer) -> bytes: ...
    def derive_into(self, key_material: Buffer, buffer: Buffer) -> int: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...
