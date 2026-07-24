from _typeshed import Incomplete, SupportsLenAndGetItem
from collections.abc import Generator, Iterable
from logging import Logger
from typing import ClassVar, Protocol, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

from .enums import AccessPermission, EncryptionMethod
from .fpdf import FPDF
from .syntax import Name, PDFObject

_Key: TypeAlias = SupportsLenAndGetItem[int]
_T_co = TypeVar("_T_co", covariant=True)

LOGGER: Logger

import_error: ImportError | None

@type_check_only
class _SupportsGetItem(Protocol[_T_co]):
    def __getitem__(self, k: int, /) -> _T_co: ...

class ARC4:
    MOD: ClassVar[int]
    def KSA(self, key: _Key) -> list[int]: ...
    def PRGA(self, S: _SupportsGetItem[int]) -> Generator[int]: ...
    def encrypt(self, key: _Key, text: Iterable[int]) -> list[int]: ...

class CryptFilter:
    type: Name
    c_f_m: Name
    length: int
    def __init__(self, mode: str, length: int) -> None: ...
    def serialize(self) -> str: ...

class EncryptionDictionary(PDFObject):
    filter: Name
    length: int
    r: int
    o: str
    u: str
    v: int
    p: int
    encrypt_metadata: str  # not always defined
    c_f: str  # not always defined
    stm_f: Name
    str_f: Name
    def __init__(self, security_handler: StandardSecurityHandler) -> None: ...

class StandardSecurityHandler:
    DEFAULT_PADDING: ClassVar[bytes]
    fpdf: FPDF
    access_permission: int
    owner_password: str
    user_password: str
    encryption_method: EncryptionMethod | None
    cf: CryptFilter | None
    key_length: int
    version: int
    revision: int
    encrypt_metadata: bool

    # The following fields are only defined after a call to generate_passwords().
    file_id: Incomplete
    info_id: Incomplete
    o: str
    k: str
    u: str
    # The following field is only defined after a call to generate_user_password_rev6().
    ue: Incomplete
    # The following field is only defined after a call to generate_owner_password_rev6().
    oe: Incomplete
    # The following field is only defined after a call to generate_perms_rev6().
    perms_rev6: Incomplete

    def __init__(
        self,
        fpdf: FPDF,
        owner_password: str,
        user_password: str | None = None,
        permission: AccessPermission = ...,
        encryption_method: EncryptionMethod = ...,
        encrypt_metadata: bool = False,
    ) -> None: ...
    def generate_passwords(self, file_id: str) -> None: ...
    def get_encryption_obj(self) -> EncryptionDictionary: ...
    @overload
    def encrypt(self, text: bytes | bytearray, obj_id: int) -> bytes: ...
    @overload
    def encrypt(self, text: str, obj_id: int) -> str: ...
    def encrypt_string(self, string: str, obj_id: int) -> str: ...
    def encrypt_stream(self, stream: bytes, obj_id: int) -> bytes: ...
    def is_aes_algorithm(self) -> bool: ...
    def encrypt_bytes(self, data: bytes, obj_id: int) -> list[int]: ...
    def encrypt_AES_cryptography(self, key: bytes, data: bytes) -> bytes: ...
    @classmethod
    def get_random_bytes(cls, size: int) -> bytes: ...
    @classmethod
    def prepare_string(cls, string: str) -> bytes: ...
    def padded_password(self, password: str) -> bytearray: ...
    def generate_owner_password(self) -> str: ...
    def generate_user_password(self) -> str: ...
    @classmethod
    def compute_hash(cls, input_password: bytes, salt: bytes, user_key: bytes = ...) -> bytes: ...
    def generate_user_password_rev6(self) -> None: ...
    def generate_owner_password_rev6(self) -> None: ...
    def generate_perms_rev6(self) -> None: ...
    def generate_encryption_key(self) -> bytes: ...

def md5(data: bytes | bytearray) -> bytes: ...
def int32(n: int) -> int: ...
