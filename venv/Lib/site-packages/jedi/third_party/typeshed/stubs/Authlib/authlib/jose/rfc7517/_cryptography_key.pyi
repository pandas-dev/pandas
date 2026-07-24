from _typeshed import ReadableBuffer
from collections.abc import Iterable
from typing import Literal, SupportsBytes, SupportsIndex, overload

from cryptography.hazmat.primitives.asymmetric.types import PrivateKeyTypes, PublicKeyTypes
from cryptography.hazmat.primitives.serialization.ssh import SSHPublicKeyTypes

@overload  # if ssh_type is None
def load_pem_key(
    raw: str | bytes | float | Iterable[SupportsIndex] | SupportsIndex | SupportsBytes | ReadableBuffer,
    ssh_type: None = None,
    key_type: str | None = None,
    password: bytes | None = None,
) -> PublicKeyTypes | PrivateKeyTypes: ...
@overload  # if key_type == "public"
def load_pem_key(
    raw: str | bytes | float | Iterable[SupportsIndex] | SupportsIndex | SupportsBytes | ReadableBuffer,
    ssh_type: ReadableBuffer | tuple[ReadableBuffer, ...] | None = None,
    key_type: Literal["public"] = ...,
    password: bytes | None = None,
) -> PublicKeyTypes: ...
@overload  # if key_type is not empty, but not "public"
def load_pem_key(
    raw: str | bytes | float | Iterable[SupportsIndex] | SupportsIndex | SupportsBytes | ReadableBuffer,
    ssh_type: ReadableBuffer | tuple[ReadableBuffer, ...] | None = None,
    key_type: str = ...,
    password: bytes | None = None,
) -> PrivateKeyTypes: ...
@overload  # if ssh_type is not empty
def load_pem_key(
    raw: str | bytes | float | Iterable[SupportsIndex] | SupportsIndex | SupportsBytes | ReadableBuffer,
    ssh_type: ReadableBuffer | tuple[ReadableBuffer, ...] = ...,
    key_type: str | None = None,
    password: bytes | None = None,
) -> SSHPublicKeyTypes | PublicKeyTypes | PrivateKeyTypes: ...
