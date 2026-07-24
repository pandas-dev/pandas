from collections.abc import Container, Iterable, Mapping, MutableMapping
from datetime import timezone
from typing import Any

from .backends.base import Key

UTC: timezone

def encode(
    claims: MutableMapping[str, Any],
    # Internally it calls jws.sign() that expects a key dict instance instead of Mapping
    key: str | bytes | dict[str, Any] | Key,
    algorithm: str = "HS256",
    headers: Mapping[str, Any] | None = None,
    access_token: str | None = None,
) -> str: ...
def decode(
    token: str | bytes,
    key: str | bytes | Mapping[str, Any] | Key | Iterable[str],
    algorithms: str | Container[str] | None = None,
    options: Mapping[str, Any] | None = None,
    audience: str | None = None,
    issuer: str | Iterable[str] | None = None,
    subject: str | None = None,
    access_token: str | None = None,
) -> dict[str, Any]: ...
def get_unverified_header(token: str | bytes) -> dict[str, Any]: ...
def get_unverified_headers(token: str | bytes) -> dict[str, Any]: ...
def get_unverified_claims(token: str | bytes) -> dict[str, Any]: ...
