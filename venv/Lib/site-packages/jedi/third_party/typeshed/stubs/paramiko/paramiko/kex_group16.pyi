from _typeshed import ReadableBuffer
from collections.abc import Callable
from hashlib import _Hash

from paramiko.kex_group1 import KexGroup1 as KexGroup1

class KexGroup16SHA512(KexGroup1):
    name: str
    P: int
    G: int
    hash_algo: Callable[[ReadableBuffer], _Hash]
