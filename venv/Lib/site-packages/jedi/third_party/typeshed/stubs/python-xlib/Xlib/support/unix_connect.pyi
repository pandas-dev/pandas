import sys
from _socket import _Address
from _typeshed import Unused
from platform import uname_result
from re import Pattern
from socket import socket
from typing import Final, Literal
from typing_extensions import TypeAlias

if sys.platform == "darwin":
    SUPPORTED_PROTOCOLS: Final[tuple[None, Literal["tcp"], Literal["unix"], Literal["darwin"]]]
    _Protocol: TypeAlias = Literal["tcp", "unix", "darwin"] | None
    DARWIN_DISPLAY_RE: Final[Pattern[str]]
else:
    SUPPORTED_PROTOCOLS: Final[tuple[None, Literal["tcp"], Literal["unix"]]]
    _Protocol: TypeAlias = Literal["tcp", "unix"] | None
uname: uname_result
DISPLAY_RE: Final[Pattern[str]]

def get_display(display: str | None) -> tuple[str, str | None, str | None, int, int]: ...
def get_socket(dname: _Address, protocol: _Protocol, host: _Address | None, dno: int) -> socket: ...
def new_get_auth(sock: socket, dname: Unused, protocol: _Protocol, host: Unused, dno: int) -> tuple[bytes, bytes]: ...
def old_get_auth(sock: Unused, dname: _Address, host: Unused, dno: Unused) -> tuple[str | Literal[b""], bytes]: ...

get_auth = new_get_auth
