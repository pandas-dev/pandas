from collections.abc import Callable, Mapping
from hashlib import _Hash
from typing import Final

from .backends.base import Key

class Algorithms:
    NONE: str
    HS256: str
    HS384: str
    HS512: str
    RS256: str
    RS384: str
    RS512: str
    ES256: str
    ES384: str
    ES512: str
    A128CBC_HS256: str
    A192CBC_HS384: str
    A256CBC_HS512: str
    A128GCM: str
    A192GCM: str
    A256GCM: str
    A128CBC: str
    A192CBC: str
    A256CBC: str
    DIR: str
    RSA1_5: str
    RSA_OAEP: str
    RSA_OAEP_256: str
    A128KW: str
    A192KW: str
    A256KW: str
    ECDH_ES: str
    ECDH_ES_A128KW: str
    ECDH_ES_A192KW: str
    ECDH_ES_A256KW: str
    A128GCMKW: str
    A192GCMKW: str
    A256GCMKW: str
    PBES2_HS256_A128KW: str
    PBES2_HS384_A192KW: str
    PBES2_HS512_A256KW: str
    DEF: str
    HMAC: set[str]
    RSA_DS: set[str]
    RSA_KW: set[str]
    RSA: set[str]
    EC_DS: set[str]
    EC_KW: set[str]
    EC: set[str]
    AES_PSEUDO: set[str]
    AES_JWE_ENC: set[str]
    AES_ENC: set[str]
    AES_KW: set[str]
    AEC_GCM_KW: set[str]
    AES: set[str]
    PBES2_KW: set[str]
    HMAC_AUTH_TAG: set[str]
    GCM: set[str]
    SUPPORTED: set[str]
    ALL: set[str]
    HASHES: Mapping[str, Callable[[bytes], _Hash]]
    KEYS: Mapping[str, type[Key]]

ALGORITHMS: Algorithms

class Zips:
    DEF: str
    NONE: None
    SUPPORTED: set[str | None]

ZIPS: Zips

JWE_SIZE_LIMIT: Final[int]
