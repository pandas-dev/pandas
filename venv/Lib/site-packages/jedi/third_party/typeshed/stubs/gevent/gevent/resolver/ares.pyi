from collections.abc import Sequence
from typing import TypedDict, type_check_only

from gevent._types import _Watcher
from gevent.hub import Hub
from gevent.resolver import AbstractResolver
from gevent.resolver.cares import channel

@type_check_only
class _ChannelArgs(TypedDict):
    flags: str | int | None
    timeout: str | float | None
    tries: str | int | None
    ndots: str | int | None
    udp_port: str | int | None
    tcp_port: str | int | None
    servers: Sequence[str] | str | None

class Resolver(AbstractResolver):
    cares_class: type[channel]
    hub: Hub
    cares: channel
    pid: int
    params: _ChannelArgs
    fork_watcher: _Watcher
    def __init__(
        self,
        hub: Hub | None = None,
        use_environ: bool = True,
        *,
        flags: str | int | None = None,
        timeout: str | float | None = None,
        tries: str | int | None = None,
        ndots: str | int | None = None,
        udp_port: str | int | None = None,
        tcp_port: str | int | None = None,
        servers: Sequence[str] | str | None = None,
    ) -> None: ...
    def __del__(self) -> None: ...

__all__ = ["Resolver"]
