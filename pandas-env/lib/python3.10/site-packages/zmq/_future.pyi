"""type annotations for async sockets"""

from __future__ import annotations

from asyncio import Future
from pickle import DEFAULT_PROTOCOL
from typing import Any, Awaitable, Literal, Sequence, TypeVar, overload

import zmq as _zmq

class _AsyncPoller(_zmq.Poller):
    _socket_class: type[_AsyncSocket]

    def poll(self, timeout=-1) -> Awaitable[list[tuple[Any, int]]]: ...  # type: ignore

T = TypeVar("T", bound="_AsyncSocket")

class _AsyncSocket(_zmq.Socket[Future]):
    @classmethod
    def from_socket(cls: type[T], socket: _zmq.Socket, io_loop: Any = None) -> T: ...
    def send(  # type: ignore
        self,
        data: Any,
        flags: int = 0,
        copy: bool = True,
        track: bool = False,
        routing_id: int | None = None,
        group: str | None = None,
    ) -> Awaitable[_zmq.MessageTracker | None]: ...
    @overload  # type: ignore
    def recv(self, flags: int = 0, *, track: bool = False) -> Awaitable[bytes]: ...
    @overload
    def recv(
        self, flags: int = 0, *, copy: Literal[True], track: bool = False
    ) -> Awaitable[bytes]: ...
    @overload
    def recv(
        self, flags: int = 0, *, copy: Literal[False], track: bool = False
    ) -> Awaitable[_zmq.Frame]: ...
    @overload
    def recv(
        self, flags: int = 0, copy: bool = True, track: bool = False
    ) -> Awaitable[bytes | _zmq.Frame]: ...
    def send_multipart(  # type: ignore
        self,
        msg_parts: Sequence,
        flags: int = 0,
        copy: bool = True,
        track: bool = False,
        routing_id: int | None = None,
        group: str | None = None,
    ) -> Awaitable[_zmq.MessageTracker | None]: ...
    @overload  # type: ignore
    def recv_multipart(
        self, flags: int = 0, *, track: bool = False
    ) -> Awaitable[list[bytes]]: ...
    @overload
    def recv_multipart(
        self, flags: int = 0, *, copy: Literal[True], track: bool = False
    ) -> Awaitable[list[bytes]]: ...
    @overload
    def recv_multipart(
        self, flags: int = 0, *, copy: Literal[False], track: bool = False
    ) -> Awaitable[list[_zmq.Frame]]: ...
    @overload
    def recv_multipart(
        self, flags: int = 0, copy: bool = True, track: bool = False
    ) -> Awaitable[list[bytes] | list[_zmq.Frame]]: ...

    # serialization wrappers

    def send_string(  # type: ignore
        self,
        u: str,
        flags: int = 0,
        copy: bool = True,
        *,
        encoding: str = 'utf-8',
        **kwargs,
    ) -> Awaitable[_zmq.Frame | None]: ...
    def recv_string(  # type: ignore
        self, flags: int = 0, encoding: str = 'utf-8'
    ) -> Awaitable[str]: ...
    def send_pyobj(  # type: ignore
        self, obj: Any, flags: int = 0, protocol: int = DEFAULT_PROTOCOL, **kwargs
    ) -> Awaitable[_zmq.Frame | None]: ...
    def recv_pyobj(self, flags: int = 0) -> Awaitable[Any]: ...  # type: ignore
    def send_json(  # type: ignore
        self, obj: Any, flags: int = 0, **kwargs
    ) -> Awaitable[_zmq.Frame | None]: ...
    def recv_json(self, flags: int = 0, **kwargs) -> Awaitable[Any]: ...  # type: ignore
    def poll(self, timeout=-1) -> Awaitable[list[tuple[Any, int]]]: ...  # type: ignore
