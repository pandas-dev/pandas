from typing import Any, Literal

from .backends import AESKey as AESKey, ECKey as ECKey, HMACKey as HMACKey, RSAKey as RSAKey
from .backends.base import DIRKey as DIRKey, Key

def get_key(algorithm: str) -> type[Key] | None: ...
def register_key(algorithm: str, key_class: type[Key]) -> Literal[True]: ...
def construct(
    # explicitly checks for key_data as dict instance, instead of a Mapping
    key_data: str | bytes | dict[str, Any] | Key,
    algorithm: str | None = None,
) -> Key: ...
