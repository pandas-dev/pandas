# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

from __future__ import annotations

from cryptography import utils
from cryptography.exceptions import UnsupportedAlgorithm, _Reasons
from cryptography.hazmat.decrepit.ciphers.modes import CFB as CFB
from cryptography.hazmat.decrepit.ciphers.modes import CFB8 as CFB8
from cryptography.hazmat.decrepit.ciphers.modes import OFB as OFB
from cryptography.hazmat.primitives._cipheralgorithm import (
    BlockCipherAlgorithm,
    CipherAlgorithm,
)
from cryptography.hazmat.primitives._modes import (
    Mode as Mode,
)
from cryptography.hazmat.primitives._modes import (
    ModeWithAuthenticationTag as ModeWithAuthenticationTag,
)
from cryptography.hazmat.primitives._modes import (
    ModeWithInitializationVector as ModeWithInitializationVector,
)
from cryptography.hazmat.primitives._modes import (
    ModeWithNonce as ModeWithNonce,
)
from cryptography.hazmat.primitives._modes import (
    ModeWithTweak as ModeWithTweak,
)
from cryptography.hazmat.primitives._modes import (
    _check_aes_key_length,
    _check_iv_and_key_length,
    _check_nonce_length,
)
from cryptography.hazmat.primitives.ciphers import algorithms


class CBC(ModeWithInitializationVector):
    name = "CBC"

    def __init__(self, initialization_vector: utils.Buffer):
        utils._check_byteslike("initialization_vector", initialization_vector)
        self._initialization_vector = initialization_vector

    @property
    def initialization_vector(self) -> utils.Buffer:
        return self._initialization_vector

    validate_for_algorithm = _check_iv_and_key_length


class XTS(ModeWithTweak):
    name = "XTS"

    def __init__(self, tweak: utils.Buffer):
        utils._check_byteslike("tweak", tweak)

        if len(tweak) != 16:
            raise ValueError("tweak must be 128-bits (16 bytes)")

        self._tweak = tweak

    @property
    def tweak(self) -> utils.Buffer:
        return self._tweak

    def validate_for_algorithm(self, algorithm: CipherAlgorithm) -> None:
        if isinstance(algorithm, (algorithms.AES128, algorithms.AES256)):
            raise TypeError(
                "The AES128 and AES256 classes do not support XTS, please use "
                "the standard AES class instead."
            )

        if algorithm.key_size not in (256, 512):
            raise ValueError(
                "The XTS specification requires a 256-bit key for AES-128-XTS"
                " and 512-bit key for AES-256-XTS"
            )


class ECB(Mode):
    name = "ECB"

    validate_for_algorithm = _check_aes_key_length


class CTR(ModeWithNonce):
    name = "CTR"

    def __init__(self, nonce: utils.Buffer):
        utils._check_byteslike("nonce", nonce)
        self._nonce = nonce

    @property
    def nonce(self) -> utils.Buffer:
        return self._nonce

    def validate_for_algorithm(self, algorithm: CipherAlgorithm) -> None:
        _check_aes_key_length(self, algorithm)
        _check_nonce_length(self.nonce, self.name, algorithm)


class GCM(ModeWithInitializationVector, ModeWithAuthenticationTag):
    name = "GCM"
    _MAX_ENCRYPTED_BYTES = (2**39 - 256) // 8
    _MAX_AAD_BYTES = (2**64) // 8

    def __init__(
        self,
        initialization_vector: utils.Buffer,
        tag: bytes | None = None,
        min_tag_length: int = 16,
    ):
        # OpenSSL 3.0.0 constrains GCM IVs to [64, 1024] bits inclusive
        # This is a sane limit anyway so we'll enforce it here.
        utils._check_byteslike("initialization_vector", initialization_vector)
        if len(initialization_vector) < 8 or len(initialization_vector) > 128:
            raise ValueError(
                "initialization_vector must be between 8 and 128 bytes (64 "
                "and 1024 bits)."
            )
        self._initialization_vector = initialization_vector
        if min_tag_length < 4:
            raise ValueError("min_tag_length must be >= 4")
        if tag is not None:
            utils._check_bytes("tag", tag)
            if len(tag) < min_tag_length:
                raise ValueError(
                    f"Authentication tag must be {min_tag_length} bytes or "
                    "longer."
                )
        self._tag = tag
        self._min_tag_length = min_tag_length

    @property
    def tag(self) -> bytes | None:
        return self._tag

    @property
    def initialization_vector(self) -> utils.Buffer:
        return self._initialization_vector

    def validate_for_algorithm(self, algorithm: CipherAlgorithm) -> None:
        _check_aes_key_length(self, algorithm)
        if not isinstance(algorithm, BlockCipherAlgorithm):
            raise UnsupportedAlgorithm(
                "GCM requires a block cipher algorithm",
                _Reasons.UNSUPPORTED_CIPHER,
            )
        block_size_bytes = algorithm.block_size // 8
        if self._tag is not None and len(self._tag) > block_size_bytes:
            raise ValueError(
                f"Authentication tag cannot be more than {block_size_bytes} "
                "bytes."
            )


utils.deprecated(
    OFB,
    __name__,
    "OFB has been moved to "
    "cryptography.hazmat.decrepit.ciphers.modes.OFB and "
    "will be removed from "
    "cryptography.hazmat.primitives.ciphers.modes in 49.0.0.",
    utils.DeprecatedIn47,
    name="OFB",
)


utils.deprecated(
    CFB,
    __name__,
    "CFB has been moved to "
    "cryptography.hazmat.decrepit.ciphers.modes.CFB and "
    "will be removed from "
    "cryptography.hazmat.primitives.ciphers.modes in 49.0.0.",
    utils.DeprecatedIn47,
    name="CFB",
)


utils.deprecated(
    CFB8,
    __name__,
    "CFB8 has been moved to "
    "cryptography.hazmat.decrepit.ciphers.modes.CFB8 and "
    "will be removed from "
    "cryptography.hazmat.primitives.ciphers.modes in 49.0.0.",
    utils.DeprecatedIn47,
    name="CFB8",
)
