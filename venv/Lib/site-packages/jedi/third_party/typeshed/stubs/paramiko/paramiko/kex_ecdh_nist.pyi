from _typeshed import ReadableBuffer
from collections.abc import Callable
from hashlib import _Hash

from cryptography.hazmat.primitives.asymmetric.ec import EllipticCurve, EllipticCurvePrivateKey, EllipticCurvePublicKey
from paramiko.message import Message
from paramiko.transport import Transport

c_MSG_KEXECDH_INIT: bytes
c_MSG_KEXECDH_REPLY: bytes

class KexNistp256:
    name: str
    hash_algo: Callable[[ReadableBuffer], _Hash]
    curve: EllipticCurve
    transport: Transport
    P: int | EllipticCurvePrivateKey
    Q_C: EllipticCurvePublicKey | None
    Q_S: EllipticCurvePublicKey | None
    def __init__(self, transport: Transport) -> None: ...
    def start_kex(self) -> None: ...
    def parse_next(self, ptype: int, m: Message) -> None: ...

class KexNistp384(KexNistp256):
    name: str
    hash_algo: Callable[[ReadableBuffer], _Hash]
    curve: EllipticCurve

class KexNistp521(KexNistp256):
    name: str
    hash_algo: Callable[[ReadableBuffer], _Hash]
    curve: EllipticCurve
