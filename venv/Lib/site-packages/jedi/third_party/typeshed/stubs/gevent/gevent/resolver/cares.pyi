from collections.abc import Callable, Iterable, Sequence
from typing import Any, Generic, TypeVar
from typing_extensions import Self, disjoint_base

from gevent._types import _AddrinfoResult, _Loop, _NameinfoResult, _SockAddr

_T = TypeVar("_T")

class ares_host_result(tuple[str, list[str], list[str]]):
    family: int
    def __new__(cls, family: int, iterable: Iterable[Any]) -> Self: ...

@disjoint_base
class Result(Generic[_T]):
    exception: BaseException | None
    value: _T | None
    def __init__(self, value: _T | None = None, exception: BaseException | None = None) -> None: ...
    def get(self) -> Any | None: ...
    def successful(self) -> bool: ...

@disjoint_base
class channel:
    @property
    def loop(self) -> _Loop: ...
    def __init__(
        self,
        loop: _Loop,
        flags: str | int | None = None,
        timeout: str | float | None = None,
        tries: str | int | None = None,
        ndots: str | int | None = None,
        udp_port: str | int | None = None,
        tcp_port: str | int | None = None,
        servers: Sequence[str] | str | None = None,
    ) -> None: ...
    def destroy(self) -> None: ...
    def getaddrinfo(
        self,
        callback: Callable[[Result[_AddrinfoResult]], object],
        name: str,
        service: str | None,
        family: int = 0,
        type: int = 0,
        proto: int = 0,
        flags: int = 0,
    ) -> None: ...
    def gethostbyaddr(self, callback: Callable[[Result[ares_host_result]], object], addr: str) -> Any: ...
    def gethostbyname(self, callback: Callable[[Result[ares_host_result]], object], name: str, family: int = 2) -> None: ...
    def getnameinfo(self, callback: Callable[[Result[_NameinfoResult]], object], sockaddr: _SockAddr, flags: int) -> None: ...
    def set_servers(self, servers: Sequence[str] | str | None = None) -> None: ...

__all__ = ["channel"]
