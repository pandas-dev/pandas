# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

import typing

from cryptography.hazmat.primitives.hashes import HashAlgorithm

def derive_pbkdf2_hmac(
    key_material: bytes,
    algorithm: HashAlgorithm,
    salt: bytes,
    iterations: int,
    length: int,
) -> bytes: ...

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
    def derive(self, key_material: bytes) -> bytes: ...
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...

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
    def verify(self, key_material: bytes, expected_key: bytes) -> None: ...
