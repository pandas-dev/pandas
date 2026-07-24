# This file is dual licensed under the terms of the Apache License, Version
# 2.0, and the BSD License. See the LICENSE file in the root of this repository
# for complete details.

from __future__ import annotations

from cryptography import utils
from cryptography.hazmat.bindings._rust import openssl as rust_openssl
from cryptography.hazmat.primitives.kdf import KeyDerivationFunction


class Mode(utils.Enum):
    CounterMode = "ctr"


class CounterLocation(utils.Enum):
    BeforeFixed = "before_fixed"
    AfterFixed = "after_fixed"
    MiddleFixed = "middle_fixed"


KBKDFHMAC = rust_openssl.kdf.KBKDFHMAC
KeyDerivationFunction.register(KBKDFHMAC)

KBKDFCMAC = rust_openssl.kdf.KBKDFCMAC
KeyDerivationFunction.register(KBKDFCMAC)
