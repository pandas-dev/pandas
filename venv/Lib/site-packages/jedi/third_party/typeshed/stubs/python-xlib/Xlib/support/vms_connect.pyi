from _socket import _Address
from _typeshed import Unused
from re import Pattern
from socket import socket
from typing import Final

display_re: Final[Pattern[str]]

def get_display(display: str | None) -> tuple[str, None, str, int, int]: ...
def get_socket(dname: _Address, protocol: Unused, host: _Address, dno: int) -> socket: ...
def get_auth(sock: Unused, dname: Unused, host: Unused, dno: Unused) -> tuple[str, str]: ...
