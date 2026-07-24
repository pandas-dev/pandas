# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

import typing

from cryptography.hazmat.primitives import padding
from cryptography.hazmat.primitives._serialization import (
    KeySerializationEncryptionBuilder,
)
from cryptography.utils import Buffer

class PKCS7PaddingContext(padding.PaddingContext):
    def __init__(self, block_size: int) -> None: ...
    def update(self, data: Buffer) -> bytes: ...
    def finalize(self) -> bytes: ...

class ANSIX923PaddingContext(padding.PaddingContext):
    def __init__(self, block_size: int) -> None: ...
    def update(self, data: Buffer) -> bytes: ...
    def finalize(self) -> bytes: ...

class PKCS7UnpaddingContext(padding.PaddingContext):
    def __init__(self, block_size: int) -> None: ...
    def update(self, data: Buffer) -> bytes: ...
    def finalize(self) -> bytes: ...

class ANSIX923UnpaddingContext(padding.PaddingContext):
    def __init__(self, block_size: int) -> None: ...
    def update(self, data: Buffer) -> bytes: ...
    def finalize(self) -> bytes: ...

class Encoding:
    PEM: typing.ClassVar[Encoding]
    DER: typing.ClassVar[Encoding]
    OpenSSH: typing.ClassVar[Encoding]
    Raw: typing.ClassVar[Encoding]
    X962: typing.ClassVar[Encoding]
    SMIME: typing.ClassVar[Encoding]

class PrivateFormat:
    PKCS8: typing.ClassVar[PrivateFormat]
    TraditionalOpenSSL: typing.ClassVar[PrivateFormat]
    Raw: typing.ClassVar[PrivateFormat]
    OpenSSH: typing.ClassVar[PrivateFormat]
    PKCS12: typing.ClassVar[PrivateFormat]
    def encryption_builder(self) -> KeySerializationEncryptionBuilder: ...

class PublicFormat:
    SubjectPublicKeyInfo: typing.ClassVar[PublicFormat]
    PKCS1: typing.ClassVar[PublicFormat]
    OpenSSH: typing.ClassVar[PublicFormat]
    Raw: typing.ClassVar[PublicFormat]
    CompressedPoint: typing.ClassVar[PublicFormat]
    UncompressedPoint: typing.ClassVar[PublicFormat]

class ParameterFormat:
    PKCS3: typing.ClassVar[ParameterFormat]

class ObjectIdentifier:
    def __init__(self, value: str) -> None: ...
    @property
    def dotted_string(self) -> str: ...
    @property
    def _name(self) -> str: ...

T = typing.TypeVar("T")
