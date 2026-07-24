# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

from __future__ import annotations

import abc

from cryptography.exceptions import UnsupportedAlgorithm, _Reasons
from cryptography.hazmat.bindings._rust import openssl as rust_openssl
from cryptography.hazmat.primitives import _serialization
from cryptography.utils import Buffer


class MLKEM768PublicKey(metaclass=abc.ABCMeta):
    @classmethod
    def from_public_bytes(cls, data: Buffer) -> MLKEM768PublicKey:
        from cryptography.hazmat.backends.openssl.backend import backend

        if not backend.mlkem_supported():
            raise UnsupportedAlgorithm(
                "ML-KEM-768 is not supported by this backend.",
                _Reasons.UNSUPPORTED_PUBLIC_KEY_ALGORITHM,
            )

        return rust_openssl.mlkem.from_mlkem768_public_bytes(data)

    @abc.abstractmethod
    def encapsulate(self) -> tuple[bytes, bytes]:
        """
        Encapsulate: returns (shared_secret, ciphertext).
        """

    @abc.abstractmethod
    def public_bytes(
        self,
        encoding: _serialization.Encoding,
        format: _serialization.PublicFormat,
    ) -> bytes:
        """
        The serialized bytes of the public key.
        """

    @abc.abstractmethod
    def public_bytes_raw(self) -> bytes:
        """
        The raw bytes of the public key.
        Equivalent to public_bytes(Raw, Raw).

        The public key is 1,184 bytes for ML-KEM-768.
        """

    @abc.abstractmethod
    def __eq__(self, other: object) -> bool:
        """
        Checks equality.
        """

    @abc.abstractmethod
    def __copy__(self) -> MLKEM768PublicKey:
        """
        Returns a copy.
        """

    @abc.abstractmethod
    def __deepcopy__(self, memo: dict) -> MLKEM768PublicKey:
        """
        Returns a deep copy.
        """


if hasattr(rust_openssl, "mlkem"):
    MLKEM768PublicKey.register(rust_openssl.mlkem.MLKEM768PublicKey)


class MLKEM768PrivateKey(metaclass=abc.ABCMeta):
    @classmethod
    def generate(cls) -> MLKEM768PrivateKey:
        from cryptography.hazmat.backends.openssl.backend import backend

        if not backend.mlkem_supported():
            raise UnsupportedAlgorithm(
                "ML-KEM-768 is not supported by this backend.",
                _Reasons.UNSUPPORTED_PUBLIC_KEY_ALGORITHM,
            )

        return rust_openssl.mlkem.generate_mlkem768_key()

    @classmethod
    def from_seed_bytes(cls, data: Buffer) -> MLKEM768PrivateKey:
        from cryptography.hazmat.backends.openssl.backend import backend

        if not backend.mlkem_supported():
            raise UnsupportedAlgorithm(
                "ML-KEM-768 is not supported by this backend.",
                _Reasons.UNSUPPORTED_PUBLIC_KEY_ALGORITHM,
            )

        return rust_openssl.mlkem.from_mlkem768_seed_bytes(data)

    @abc.abstractmethod
    def decapsulate(self, ciphertext: Buffer) -> bytes:
        """
        Decapsulate: returns shared_secret.
        """

    @abc.abstractmethod
    def public_key(self) -> MLKEM768PublicKey:
        """
        The MLKEM768PublicKey derived from this private key.
        """

    @abc.abstractmethod
    def private_bytes(
        self,
        encoding: _serialization.Encoding,
        format: _serialization.PrivateFormat,
        encryption_algorithm: _serialization.KeySerializationEncryption,
    ) -> bytes:
        """
        The serialized bytes of the private key.
        """

    @abc.abstractmethod
    def private_bytes_raw(self) -> bytes:
        """
        The raw bytes of the private key (64-byte seed).
        Equivalent to private_bytes(Raw, Raw, NoEncryption()).
        """

    @abc.abstractmethod
    def __copy__(self) -> MLKEM768PrivateKey:
        """
        Returns a copy.
        """

    @abc.abstractmethod
    def __deepcopy__(self, memo: dict) -> MLKEM768PrivateKey:
        """
        Returns a deep copy.
        """


if hasattr(rust_openssl, "mlkem"):
    MLKEM768PrivateKey.register(rust_openssl.mlkem.MLKEM768PrivateKey)


class MLKEM1024PublicKey(metaclass=abc.ABCMeta):
    @classmethod
    def from_public_bytes(cls, data: Buffer) -> MLKEM1024PublicKey:
        from cryptography.hazmat.backends.openssl.backend import backend

        if not backend.mlkem_supported():
            raise UnsupportedAlgorithm(
                "ML-KEM-1024 is not supported by this backend.",
                _Reasons.UNSUPPORTED_PUBLIC_KEY_ALGORITHM,
            )

        return rust_openssl.mlkem.from_mlkem1024_public_bytes(data)

    @abc.abstractmethod
    def encapsulate(self) -> tuple[bytes, bytes]:
        """
        Encapsulate: returns (shared_secret, ciphertext).
        """

    @abc.abstractmethod
    def public_bytes(
        self,
        encoding: _serialization.Encoding,
        format: _serialization.PublicFormat,
    ) -> bytes:
        """
        The serialized bytes of the public key.
        """

    @abc.abstractmethod
    def public_bytes_raw(self) -> bytes:
        """
        The raw bytes of the public key.
        Equivalent to public_bytes(Raw, Raw).

        The public key is 1,568 bytes for ML-KEM-1024.
        """

    @abc.abstractmethod
    def __eq__(self, other: object) -> bool:
        """
        Checks equality.
        """

    @abc.abstractmethod
    def __copy__(self) -> MLKEM1024PublicKey:
        """
        Returns a copy.
        """

    @abc.abstractmethod
    def __deepcopy__(self, memo: dict) -> MLKEM1024PublicKey:
        """
        Returns a deep copy.
        """


if hasattr(rust_openssl, "mlkem"):
    MLKEM1024PublicKey.register(rust_openssl.mlkem.MLKEM1024PublicKey)


class MLKEM1024PrivateKey(metaclass=abc.ABCMeta):
    @classmethod
    def generate(cls) -> MLKEM1024PrivateKey:
        from cryptography.hazmat.backends.openssl.backend import backend

        if not backend.mlkem_supported():
            raise UnsupportedAlgorithm(
                "ML-KEM-1024 is not supported by this backend.",
                _Reasons.UNSUPPORTED_PUBLIC_KEY_ALGORITHM,
            )

        return rust_openssl.mlkem.generate_mlkem1024_key()

    @classmethod
    def from_seed_bytes(cls, data: Buffer) -> MLKEM1024PrivateKey:
        from cryptography.hazmat.backends.openssl.backend import backend

        if not backend.mlkem_supported():
            raise UnsupportedAlgorithm(
                "ML-KEM-1024 is not supported by this backend.",
                _Reasons.UNSUPPORTED_PUBLIC_KEY_ALGORITHM,
            )

        return rust_openssl.mlkem.from_mlkem1024_seed_bytes(data)

    @abc.abstractmethod
    def decapsulate(self, ciphertext: Buffer) -> bytes:
        """
        Decapsulate: returns shared_secret.
        """

    @abc.abstractmethod
    def public_key(self) -> MLKEM1024PublicKey:
        """
        The MLKEM1024PublicKey derived from this private key.
        """

    @abc.abstractmethod
    def private_bytes(
        self,
        encoding: _serialization.Encoding,
        format: _serialization.PrivateFormat,
        encryption_algorithm: _serialization.KeySerializationEncryption,
    ) -> bytes:
        """
        The serialized bytes of the private key.
        """

    @abc.abstractmethod
    def private_bytes_raw(self) -> bytes:
        """
        The raw bytes of the private key (64-byte seed).
        Equivalent to private_bytes(Raw, Raw, NoEncryption()).
        """

    @abc.abstractmethod
    def __copy__(self) -> MLKEM1024PrivateKey:
        """
        Returns a copy.
        """

    @abc.abstractmethod
    def __deepcopy__(self, memo: dict) -> MLKEM1024PrivateKey:
        """
        Returns a deep copy.
        """


if hasattr(rust_openssl, "mlkem"):
    MLKEM1024PrivateKey.register(rust_openssl.mlkem.MLKEM1024PrivateKey)
