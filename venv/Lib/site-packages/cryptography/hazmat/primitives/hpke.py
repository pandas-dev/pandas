# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

from __future__ import annotations

from cryptography.hazmat.bindings._rust import openssl as rust_openssl

AEAD = rust_openssl.hpke.AEAD
KDF = rust_openssl.hpke.KDF
KEM = rust_openssl.hpke.KEM
MLKEM768X25519PrivateKey = rust_openssl.hpke.MLKEM768X25519PrivateKey
MLKEM768X25519PublicKey = rust_openssl.hpke.MLKEM768X25519PublicKey
MLKEM1024P384PrivateKey = rust_openssl.hpke.MLKEM1024P384PrivateKey
MLKEM1024P384PublicKey = rust_openssl.hpke.MLKEM1024P384PublicKey
Suite = rust_openssl.hpke.Suite

__all__ = [
    "AEAD",
    "KDF",
    "KEM",
    "MLKEM768X25519PrivateKey",
    "MLKEM768X25519PublicKey",
    "MLKEM1024P384PrivateKey",
    "MLKEM1024P384PublicKey",
    "Suite",
]
