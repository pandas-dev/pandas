from _socket import SocketType
from logging import Logger

LOGGER: Logger

_SUPPORTED_TCP_OPTIONS: dict[str, int]

def socket_requires_keepalive(tcp_options: dict[str, int]) -> bool: ...
def set_sock_opts(tcp_options: dict[str, int] | None, sock: SocketType) -> None: ...
