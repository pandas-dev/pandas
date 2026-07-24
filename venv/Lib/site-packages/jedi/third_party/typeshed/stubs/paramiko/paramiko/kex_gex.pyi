from _typeshed import ReadableBuffer
from collections.abc import Callable
from hashlib import _Hash

from paramiko.message import Message
from paramiko.transport import Transport

c_MSG_KEXDH_GEX_REQUEST_OLD: bytes
c_MSG_KEXDH_GEX_GROUP: bytes
c_MSG_KEXDH_GEX_INIT: bytes
c_MSG_KEXDH_GEX_REPLY: bytes
c_MSG_KEXDH_GEX_REQUEST: bytes

class KexGex:
    name: str
    min_bits: int
    max_bits: int
    preferred_bits: int
    hash_algo: Callable[[ReadableBuffer], _Hash]
    transport: Transport
    p: int | None
    q: int | None
    g: int | None
    x: int | None
    e: int | None
    f: int | None
    old_style: bool
    def __init__(self, transport: Transport) -> None: ...
    def start_kex(self, _test_old_style: bool = False) -> None: ...
    def parse_next(self, ptype: int, m: Message) -> None: ...

class KexGexSHA256(KexGex):
    name: str
    hash_algo: Callable[[ReadableBuffer], _Hash]
