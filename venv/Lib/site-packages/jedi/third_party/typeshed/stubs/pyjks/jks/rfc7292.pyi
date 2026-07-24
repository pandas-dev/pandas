from hashlib import _Hash
from typing import Final, Literal
from typing_extensions import TypeAlias

from pyasn1.type.namedtype import NamedTypes
from pyasn1.type.univ import Sequence

PBE_WITH_SHA1_AND_TRIPLE_DES_CBC_OID: Final[tuple[int, ...]]
PURPOSE_KEY_MATERIAL: Final = 1
PURPOSE_IV_MATERIAL: Final = 2
PURPOSE_MAC_MATERIAL: Final = 3

_Purpose: TypeAlias = Literal[1, 2, 3]

class Pkcs12PBEParams(Sequence):
    componentType: NamedTypes

def derive_key(
    hashfn: _Hash, purpose_byte: _Purpose, password_str: str, salt: bytes, iteration_count: int, desired_key_size: int
) -> bytes: ...
def decrypt_PBEWithSHAAnd3KeyTripleDESCBC(
    data: bytes | bytearray, password_str: str, salt: bytes, iteration_count: int
) -> bytes: ...
def decrypt_PBEWithSHAAndTwofishCBC(
    encrypted_data: bytes | bytearray, password: str, salt: bytes, iteration_count: int
) -> bytes: ...
def encrypt_PBEWithSHAAndTwofishCBC(
    plaintext_data: bytes | bytearray, password: str, salt: bytes, iteration_count: int
) -> bytes: ...
