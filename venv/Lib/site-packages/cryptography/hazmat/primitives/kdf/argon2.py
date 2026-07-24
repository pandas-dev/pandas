# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

from __future__ import annotations

from cryptography.hazmat.bindings._rust import openssl as rust_openssl
from cryptography.hazmat.primitives.kdf import KeyDerivationFunction

Argon2d = rust_openssl.kdf.Argon2d
Argon2i = rust_openssl.kdf.Argon2i
Argon2id = rust_openssl.kdf.Argon2id
KeyDerivationFunction.register(Argon2d)
KeyDerivationFunction.register(Argon2i)
KeyDerivationFunction.register(Argon2id)

__all__ = ["Argon2d", "Argon2i", "Argon2id"]
