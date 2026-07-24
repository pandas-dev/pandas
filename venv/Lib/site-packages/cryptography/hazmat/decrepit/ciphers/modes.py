# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

from __future__ import annotations

from cryptography import utils
from cryptography.hazmat.primitives._modes import (
    ModeWithInitializationVector,
    _check_iv_and_key_length,
)


class OFB(ModeWithInitializationVector):
    name = "OFB"

    def __init__(self, initialization_vector: utils.Buffer):
        utils._check_byteslike("initialization_vector", initialization_vector)
        self._initialization_vector = initialization_vector

    @property
    def initialization_vector(self) -> utils.Buffer:
        return self._initialization_vector

    validate_for_algorithm = _check_iv_and_key_length


class CFB(ModeWithInitializationVector):
    name = "CFB"

    def __init__(self, initialization_vector: utils.Buffer):
        utils._check_byteslike("initialization_vector", initialization_vector)
        self._initialization_vector = initialization_vector

    @property
    def initialization_vector(self) -> utils.Buffer:
        return self._initialization_vector

    validate_for_algorithm = _check_iv_and_key_length


class CFB8(ModeWithInitializationVector):
    name = "CFB8"

    def __init__(self, initialization_vector: utils.Buffer):
        utils._check_byteslike("initialization_vector", initialization_vector)
        self._initialization_vector = initialization_vector

    @property
    def initialization_vector(self) -> utils.Buffer:
        return self._initialization_vector

    validate_for_algorithm = _check_iv_and_key_length
