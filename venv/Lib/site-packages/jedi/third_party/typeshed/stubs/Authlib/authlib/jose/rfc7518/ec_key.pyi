from typing import ClassVar

from authlib.jose.rfc7517 import AsymmetricKey
from cryptography.hazmat.primitives.asymmetric.ec import EllipticCurvePrivateKey, EllipticCurvePublicKey

class ECKey(AsymmetricKey):
    kty: str
    DSS_CURVES: dict[str, type]
    CURVES_DSS: dict[property, str]
    REQUIRED_JSON_FIELDS: ClassVar[list[str]]
    PUBLIC_KEY_FIELDS = REQUIRED_JSON_FIELDS
    PRIVATE_KEY_FIELDS: ClassVar[list[str]]
    PUBLIC_KEY_CLS: ClassVar[type]
    PRIVATE_KEY_CLS: ClassVar[type]
    SSH_PUBLIC_PREFIX: ClassVar[bytes]
    def exchange_shared_key(self, pubkey): ...
    @property
    def curve_key_size(self): ...
    def load_private_key(self) -> EllipticCurvePrivateKey: ...
    def load_public_key(self) -> EllipticCurvePublicKey: ...
    def dumps_private_key(self) -> dict[str, str]: ...
    def dumps_public_key(self) -> dict[str, str]: ...
    @classmethod
    def generate_key(cls, crv: str = "P-256", options=None, is_private: bool = False) -> ECKey: ...
