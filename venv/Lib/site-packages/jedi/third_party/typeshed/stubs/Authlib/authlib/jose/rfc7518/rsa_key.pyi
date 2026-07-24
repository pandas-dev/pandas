from collections.abc import Iterable
from typing import ClassVar
from typing_extensions import Self

from authlib.jose.rfc7517 import AsymmetricKey
from cryptography.hazmat.primitives.asymmetric.rsa import RSAPrivateKey, RSAPublicKey

class RSAKey(AsymmetricKey):
    kty: str
    PUBLIC_KEY_CLS: ClassVar[type]
    PRIVATE_KEY_CLS: ClassVar[type]
    PUBLIC_KEY_FIELDS: ClassVar[list[str]]
    PRIVATE_KEY_FIELDS: ClassVar[list[str]]
    REQUIRED_JSON_FIELDS: ClassVar[list[str]]
    SSH_PUBLIC_PREFIX: ClassVar[bytes]
    def dumps_private_key(self) -> dict[str, str]: ...
    def dumps_public_key(self) -> dict[str, str]: ...
    def load_private_key(self) -> RSAPrivateKey: ...
    def load_public_key(self) -> RSAPublicKey: ...
    @classmethod
    def generate_key(cls, key_size: int = 2048, options=None, is_private: bool = False) -> RSAKey: ...
    @classmethod
    def import_dict_key(cls, raw, options=None) -> Self: ...

def has_all_prime_factors(obj: Iterable[str]) -> bool: ...
