"""Types for aiohappyeyeballs."""

import socket
from collections.abc import Callable

AddrInfoType = tuple[
    int | socket.AddressFamily,
    int | socket.SocketKind,
    int,
    str,
    tuple,  # type: ignore[type-arg]
]

SocketFactoryType = Callable[[AddrInfoType], socket.socket]
