from _typeshed import Incomplete
from typing import Final

from authlib.jose.rfc7516 import JWEEncAlgorithm

class CBCHS2EncAlgorithm(JWEEncAlgorithm):
    IV_SIZE: int
    name: str
    description: str
    key_size: int
    key_len: int
    CEK_SIZE: int
    hash_alg: Incomplete
    def __init__(self, key_size: int, hash_type: int | str) -> None: ...
    def encrypt(self, msg, aad, iv, key) -> tuple[bytes, bytes]: ...
    def decrypt(self, ciphertext, aad, iv, tag, key) -> bytes: ...

class GCMEncAlgorithm(JWEEncAlgorithm):
    IV_SIZE: int
    name: str
    description: str
    key_size: int
    CEK_SIZE: int
    def __init__(self, key_size: int) -> None: ...
    def encrypt(self, msg, aad, iv, key) -> tuple[bytes, bytes]: ...
    def decrypt(self, ciphertext, aad, iv, tag, key) -> bytes: ...

JWE_ENC_ALGORITHMS: Final[list[CBCHS2EncAlgorithm | GCMEncAlgorithm]]
