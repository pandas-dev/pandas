# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

from __future__ import annotations

from cryptography.hazmat.bindings._rust import openssl as rust_openssl
from cryptography.hazmat.primitives.kdf import KeyDerivationFunction

ConcatKDFHash = rust_openssl.kdf.ConcatKDFHash
ConcatKDFHMAC = rust_openssl.kdf.ConcatKDFHMAC

KeyDerivationFunction.register(ConcatKDFHash)
KeyDerivationFunction.register(ConcatKDFHMAC)

__all__ = ["ConcatKDFHMAC", "ConcatKDFHash"]
