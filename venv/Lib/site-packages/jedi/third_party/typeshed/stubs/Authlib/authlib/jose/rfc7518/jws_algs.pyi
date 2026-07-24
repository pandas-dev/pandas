import hashlib
from _typeshed import Incomplete

from authlib.jose.rfc7515 import JWSAlgorithm
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.asymmetric.padding import PKCS1v15

from .ec_key import ECKey
from .oct_key import OctKey
from .rsa_key import RSAKey

class NoneAlgorithm(JWSAlgorithm):
    name: str
    description: str
    deprecated: bool
    def prepare_key(self, raw_data) -> None: ...
    def sign(self, msg, key) -> bytes: ...
    def verify(self, msg, sig, key) -> bool: ...

class HMACAlgorithm(JWSAlgorithm):
    SHA256 = hashlib.sha256
    SHA384 = hashlib.sha384
    SHA512 = hashlib.sha512
    name: str
    description: str
    hash_alg: Incomplete
    def __init__(self, sha_type: int | str) -> None: ...
    def prepare_key(self, raw_data) -> OctKey: ...
    def sign(self, msg, key) -> bytes: ...
    def verify(self, msg, sig, key) -> bool: ...

class RSAAlgorithm(JWSAlgorithm):
    SHA256 = hashes.SHA256
    SHA384 = hashes.SHA384
    SHA512 = hashes.SHA512
    name: str
    description: str
    hash_alg: Incomplete
    padding: PKCS1v15
    def __init__(self, sha_type: int | str) -> None: ...
    def prepare_key(self, raw_data) -> RSAKey: ...
    def sign(self, msg, key): ...
    def verify(self, msg, sig, key) -> bool: ...

class ECAlgorithm(JWSAlgorithm):
    SHA256 = hashes.SHA256
    SHA384 = hashes.SHA384
    SHA512 = hashes.SHA512
    name: str
    curve: Incomplete
    description: str
    hash_alg: Incomplete
    def __init__(self, name: str, curve, sha_type: int | str) -> None: ...
    def prepare_key(self, raw_data) -> ECKey: ...
    def sign(self, msg, key) -> bytes: ...
    def verify(self, msg, sig, key) -> bool: ...

class RSAPSSAlgorithm(JWSAlgorithm):
    SHA256 = hashes.SHA256
    SHA384 = hashes.SHA384
    SHA512 = hashes.SHA512
    name: str
    description: str
    hash_alg: Incomplete
    def __init__(self, sha_type: int | str) -> None: ...
    def prepare_key(self, raw_data) -> RSAKey: ...
    def sign(self, msg, key): ...
    def verify(self, msg, sig, key) -> bool: ...

JWS_ALGORITHMS: list[JWSAlgorithm]
