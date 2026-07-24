# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

from __future__ import annotations

import abc

from cryptography import utils
from cryptography.exceptions import UnsupportedAlgorithm, _Reasons
from cryptography.hazmat.primitives._cipheralgorithm import (
    BlockCipherAlgorithm,
    CipherAlgorithm,
)


class Mode(metaclass=abc.ABCMeta):
    @property
    @abc.abstractmethod
    def name(self) -> str:
        """
        A string naming this mode (e.g. "ECB", "CBC").
        """

    @abc.abstractmethod
    def validate_for_algorithm(self, algorithm: CipherAlgorithm) -> None:
        """
        Checks that all the necessary invariants of this (mode, algorithm)
        combination are met.
        """


class ModeWithInitializationVector(Mode, metaclass=abc.ABCMeta):
    @property
    @abc.abstractmethod
    def initialization_vector(self) -> utils.Buffer:
        """
        The value of the initialization vector for this mode as bytes.
        """


class ModeWithTweak(Mode, metaclass=abc.ABCMeta):
    @property
    @abc.abstractmethod
    def tweak(self) -> utils.Buffer:
        """
        The value of the tweak for this mode as bytes.
        """


class ModeWithNonce(Mode, metaclass=abc.ABCMeta):
    @property
    @abc.abstractmethod
    def nonce(self) -> utils.Buffer:
        """
        The value of the nonce for this mode as bytes.
        """


class ModeWithAuthenticationTag(Mode, metaclass=abc.ABCMeta):
    @property
    @abc.abstractmethod
    def tag(self) -> bytes | None:
        """
        The value of the tag supplied to the constructor of this mode.
        """


def _check_aes_key_length(self: Mode, algorithm: CipherAlgorithm) -> None:
    if algorithm.key_size > 256 and algorithm.name == "AES":
        raise ValueError(
            "Only 128, 192, and 256 bit keys are allowed for this AES mode"
        )


def _check_iv_length(
    self: ModeWithInitializationVector, algorithm: BlockCipherAlgorithm
) -> None:
    iv_len = len(self.initialization_vector)
    if iv_len * 8 != algorithm.block_size:
        raise ValueError(f"Invalid IV size ({iv_len}) for {self.name}.")


def _check_nonce_length(
    nonce: utils.Buffer, name: str, algorithm: CipherAlgorithm
) -> None:
    if not isinstance(algorithm, BlockCipherAlgorithm):
        raise UnsupportedAlgorithm(
            f"{name} requires a block cipher algorithm",
            _Reasons.UNSUPPORTED_CIPHER,
        )
    if len(nonce) * 8 != algorithm.block_size:
        raise ValueError(f"Invalid nonce size ({len(nonce)}) for {name}.")


def _check_iv_and_key_length(
    self: ModeWithInitializationVector, algorithm: CipherAlgorithm
) -> None:
    if not isinstance(algorithm, BlockCipherAlgorithm):
        raise UnsupportedAlgorithm(
            f"{self} requires a block cipher algorithm",
            _Reasons.UNSUPPORTED_CIPHER,
        )
    _check_aes_key_length(self, algorithm)
    _check_iv_length(self, algorithm)
