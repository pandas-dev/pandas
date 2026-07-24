from _typeshed import ReadableBuffer
from collections.abc import Callable
from hashlib import _Hash

from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey
from paramiko.message import Message
from paramiko.transport import Transport

c_MSG_KEXECDH_INIT: bytes
c_MSG_KEXECDH_REPLY: bytes

class KexCurve25519:
    hash_algo: Callable[[ReadableBuffer], _Hash]
    transport: Transport
    key: X25519PrivateKey | None
    def __init__(self, transport: Transport) -> None: ...
    @classmethod
    def is_available(cls) -> bool: ...
    def start_kex(self) -> None: ...
    def parse_next(self, ptype: int, m: Message) -> None: ...
