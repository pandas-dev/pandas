from _typeshed import SupportsKeysAndGetItem, Unused
from typing import Final, Literal
from typing_extensions import Self, TypeAlias

from .jks import TrustedCertEntry
from .util import AbstractKeystore, AbstractKeystoreEntry

_BksType: TypeAlias = Literal["bks", "uber"]
_CertType: TypeAlias = Literal["X.509"]
_EntryFormat: TypeAlias = Literal["PKCS8", "PKCS#8", "X.509", "X509", "RAW"]
_BksVersion: TypeAlias = Literal[1, 2]

ENTRY_TYPE_CERTIFICATE: Final = 1
ENTRY_TYPE_KEY: Final = 2
ENTRY_TYPE_SECRET: Final = 3
ENTRY_TYPE_SEALED: Final = 4

KEY_TYPE_PRIVATE: Final = 0
KEY_TYPE_PUBLIC: Final = 1
KEY_TYPE_SECRET: Final = 2
_KeyType: TypeAlias = Literal[0, 1, 2]

class AbstractBksEntry(AbstractKeystoreEntry):
    store_type: _BksType | None
    cert_chain: list[tuple[_CertType, bytes]]
    def __init__(
        self,
        *,
        cert_chain: list[tuple[_CertType, bytes]] = ...,
        encrypted: bytes | None = None,
        store_type: _BksType | None = None,
        alias: str,
        timestamp: int,
        **kwargs: Unused,
    ) -> None: ...

class BksTrustedCertEntry(TrustedCertEntry):
    store_type: _BksType | None  # type: ignore[assignment]

class BksKeyEntry(AbstractBksEntry):
    type: _KeyType
    format: _EntryFormat
    algorithm: str
    encoded: bytes
    # type == KEY_TYPE_PRIVATE
    pkey_pkcs8: bytes
    pkey: bytes
    algorithm_oid: tuple[int, ...]
    # type == KEY_TYPE_PUBLIC
    public_key_info: bytes
    public_key: bytes
    # type == KEY_TYPE_SECRET
    key: bytes
    key_size: int
    def __init__(
        self,
        type: _KeyType,
        format: _EntryFormat,
        algorithm: str,
        encoded: bytes,
        *,
        cert_chain: list[tuple[_CertType, bytes]] = ...,
        encrypted: bytes | None = None,
        store_type: _BksType | None = None,
        alias: str,
        timestamp: int,
        **kwargs: Unused,
    ) -> None: ...
    @classmethod
    def type2str(cls, t: _KeyType) -> Literal["PRIVATE", "PUBLIC", "SECRET"]: ...
    def is_decrypted(self) -> Literal[True]: ...

class BksSecretKeyEntry(AbstractBksEntry):
    key: bytes
    def is_decrypted(self) -> Literal[True]: ...

class BksSealedKeyEntry(AbstractBksEntry):
    # Properties provided by __getattr__
    nested: BksKeyEntry | None
    # __getattr__ proxies all attributes of nested BksKeyEntry after decrypting
    type: _KeyType
    format: _EntryFormat
    algorithm: str
    encoded: bytes
    # if type == KEY_TYPE_PRIVATE
    pkey_pkcs8: bytes
    pkey: bytes
    algorithm_oid: tuple[int, ...]
    # if type == KEY_TYPE_PUBLIC
    public_key_info: bytes
    public_key: bytes
    # if type == KEY_TYPE_SECRET
    key: bytes
    key_size: int

class BksKeyStore(AbstractKeystore):
    store_type: Literal["bks"]
    entries: dict[str, BksTrustedCertEntry | BksKeyEntry | BksSealedKeyEntry | BksSecretKeyEntry]  # type: ignore[assignment]
    version: _BksVersion
    def __init__(
        self,
        store_type: Literal["bks"],
        entries: SupportsKeysAndGetItem[str, BksTrustedCertEntry | BksKeyEntry | BksSealedKeyEntry | BksSecretKeyEntry],
        version: _BksVersion = 2,
    ) -> None: ...
    @property
    def certs(self) -> dict[str, BksTrustedCertEntry]: ...
    @property
    def plain_keys(self) -> dict[str, BksKeyEntry]: ...
    @property
    def sealed_keys(self) -> dict[str, BksSealedKeyEntry]: ...
    @property
    def secret_keys(self) -> dict[str, BksSecretKeyEntry]: ...
    @classmethod
    def loads(cls, data: bytes, store_password: str, try_decrypt_keys: bool = True) -> Self: ...

class UberKeyStore(BksKeyStore):
    store_type: Literal["uber"]  # type: ignore[assignment]
    version: Literal[1]
    def __init__(
        self,
        store_type: Literal["uber"],
        entries: SupportsKeysAndGetItem[str, BksTrustedCertEntry | BksKeyEntry | BksSealedKeyEntry | BksSecretKeyEntry],
        version: Literal[1] = 1,
    ) -> None: ...
