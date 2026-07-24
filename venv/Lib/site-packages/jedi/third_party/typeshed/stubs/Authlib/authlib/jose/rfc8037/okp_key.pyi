from typing import ClassVar, Final

from authlib.jose.rfc7517 import AsymmetricKey
from cryptography.hazmat.primitives.asymmetric.ed448 import Ed448PrivateKey, Ed448PublicKey
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey, Ed25519PublicKey
from cryptography.hazmat.primitives.asymmetric.x448 import X448PrivateKey, X448PublicKey
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey, X25519PublicKey

PUBLIC_KEYS_MAP: Final[dict[str, type]]
PRIVATE_KEYS_MAP: Final[dict[str, type]]

class OKPKey(AsymmetricKey):
    kty: str
    REQUIRED_JSON_FIELDS: ClassVar[list[str]]
    PUBLIC_KEY_FIELDS = REQUIRED_JSON_FIELDS
    PRIVATE_KEY_FIELDS: ClassVar[list[str]]
    PUBLIC_KEY_CLS: ClassVar[tuple[type, ...]]
    PRIVATE_KEY_CLS: ClassVar[tuple[type, ...]]
    SSH_PUBLIC_PREFIX: ClassVar[bytes]
    def exchange_shared_key(self, pubkey: X25519PrivateKey | X448PublicKey) -> bytes: ...
    @staticmethod
    def get_key_curve(key) -> str | None: ...
    def load_private_key(self) -> Ed25519PrivateKey | Ed448PrivateKey | X25519PrivateKey | X448PrivateKey: ...
    def load_public_key(self) -> Ed25519PublicKey | Ed448PublicKey | X25519PublicKey | X448PublicKey: ...
    def dumps_private_key(self) -> dict[str, str | None]: ...
    def dumps_public_key(self, public_key=None) -> dict[str, str | None]: ...
    @classmethod
    def generate_key(cls, crv: str = "Ed25519", options=None, is_private: bool = False) -> OKPKey: ...
