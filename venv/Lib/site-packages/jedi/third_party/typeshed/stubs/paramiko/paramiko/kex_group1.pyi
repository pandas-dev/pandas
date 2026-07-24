from _typeshed import ReadableBuffer
from collections.abc import Callable
from hashlib import _Hash

from paramiko.message import Message
from paramiko.transport import Transport

c_MSG_KEXDH_INIT: bytes
c_MSG_KEXDH_REPLY: bytes
b7fffffffffffffff: bytes
b0000000000000000: bytes

class KexGroup1:
    P: int
    G: int
    name: str
    hash_algo: Callable[[ReadableBuffer], _Hash]
    transport: Transport
    x: int
    e: int
    f: int
    def __init__(self, transport: Transport) -> None: ...
    def start_kex(self) -> None: ...
    def parse_next(self, ptype: int, m: Message) -> None: ...
