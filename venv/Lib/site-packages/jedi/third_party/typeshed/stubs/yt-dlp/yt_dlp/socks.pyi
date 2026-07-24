import socket
from _socket import _Address
from collections.abc import Mapping
from typing import Final, NamedTuple, SupportsIndex

SOCKS4_VERSION: Final = 4
SOCKS4_REPLY_VERSION: Final = 0x00
SOCKS4_DEFAULT_DSTIP: Final[bytes]
SOCKS5_VERSION: Final = 5
SOCKS5_USER_AUTH_VERSION: Final = 0x01
SOCKS5_USER_AUTH_SUCCESS: Final = 0x00

class Socks4Command:
    CMD_CONNECT: int
    CMD_BIND: int

class Socks5Command(Socks4Command):
    CMD_UDP_ASSOCIATE: int

class Socks5Auth:
    AUTH_NONE: int
    AUTH_GSSAPI: int
    AUTH_USER_PASS: int
    AUTH_NO_ACCEPTABLE: int

class Socks5AddressType:
    ATYP_IPV4: int
    ATYP_DOMAINNAME: int
    ATYP_IPV6: int

class ProxyError(OSError):
    ERR_SUCCESS: int
    def __init__(self, code: int | None = None, msg: str | None = None) -> None: ...

class InvalidVersionError(ProxyError):
    def __init__(self, expected_version: int, got_version: int) -> None: ...

class Socks4Error(ProxyError):
    ERR_SUCCESS: int
    CODES: Mapping[int, str]

class Socks5Error(ProxyError):
    ERR_GENERAL_FAILURE: int
    CODES: Mapping[int, str]

class ProxyType:
    SOCKS4: int
    SOCKS4A: int
    SOCKS5: int

class Proxy(NamedTuple):
    type: ProxyType
    host: str
    port: int
    username: str
    password: str
    remote_dns: bool

class sockssocket(socket.socket):
    def __init__(
        self, family: int = -1, type: int = -1, proto: int = -1, fileno: SupportsIndex | bytes | None = None
    ) -> None: ...
    def setproxy(
        self,
        proxytype: ProxyType,
        addr: str,
        port: int,
        rdns: bool = True,
        username: str | None = None,
        password: str | None = None,
    ) -> None: ...
    def recvall(self, cnt: int) -> bytes: ...
    def connect(self, address: _Address) -> None: ...
    def connect_ex(self, address: _Address) -> int: ...
