from _typeshed import Incomplete
from typing import ClassVar

class SignatureVerifier:
    DISABLE_JWT_CHECKS: ClassVar[dict[str, bool]]
    def __init__(self, algorithm: str) -> None: ...
    async def verify_signature(self, token: str) -> dict[str, Incomplete]: ...

class SymmetricSignatureVerifier(SignatureVerifier):
    def __init__(self, shared_secret: str, algorithm: str = "HS256") -> None: ...

class JwksFetcher:
    CACHE_TTL: ClassVar[int]
    def __init__(self, jwks_url: str, cache_ttl: int = 600) -> None: ...
    def get_key(self, key_id: str): ...

class AsymmetricSignatureVerifier(SignatureVerifier):
    def __init__(self, jwks_url: str, algorithm: str = "RS256", cache_ttl: int = 600) -> None: ...

class TokenVerifier:
    iss: str
    aud: str
    leeway: int
    def __init__(self, signature_verifier: SignatureVerifier, issuer: str, audience: str, leeway: int = 0) -> None: ...
    def verify(
        self, token: str, nonce: str | None = None, max_age: int | None = None, organization: str | None = None
    ) -> dict[str, Incomplete]: ...
