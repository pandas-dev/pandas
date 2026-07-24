from _typeshed import Incomplete
from collections.abc import Collection, Mapping

from authlib.jose.rfc7517 import Key, KeySet

class JsonWebKey:
    JWK_KEY_CLS: dict[Incomplete, Incomplete]
    @classmethod
    def generate_key(cls, kty, crv_or_size, options=None, is_private: bool = False): ...
    @classmethod
    def import_key(cls, raw: Mapping[str, object], options: Mapping[str, object] | None = None) -> Key: ...
    @classmethod
    def import_key_set(cls, raw: str | Collection[str] | dict[str, object]) -> KeySet: ...
