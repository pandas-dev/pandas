from _typeshed import ReadableBuffer
from collections.abc import Callable
from hashlib import _Hash

from paramiko.kex_group1 import KexGroup1 as KexGroup1

class KexGroup14(KexGroup1):
    P: int
    G: int
    name: str
    hash_algo: Callable[[ReadableBuffer], _Hash]

class KexGroup14SHA256(KexGroup14):
    name: str
    hash_algo: Callable[[ReadableBuffer], _Hash]
