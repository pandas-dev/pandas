from _typeshed import Incomplete
from collections.abc import Callable
from logging import Logger
from threading import Thread
from typing import Any, Literal, Protocol, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import ParamSpec, TypeAlias, Unpack

from flask import Flask
from flask.testing import FlaskClient

from .namespace import Namespace as Namespace
from .test_client import SocketIOTestClient as SocketIOTestClient

_P = ParamSpec("_P")
_R_co = TypeVar("_R_co", covariant=True)
_ExceptionHandler: TypeAlias = Callable[[BaseException], _R_co]
_Handler: TypeAlias = Callable[_P, _R_co]

@type_check_only
class _HandlerDecorator(Protocol):
    def __call__(self, handler: _Handler[_P, _R_co]) -> _Handler[_P, _R_co]: ...

@type_check_only
class _ExceptionHandlerDecorator(Protocol):
    def __call__(self, exception_handler: _ExceptionHandler[_R_co]) -> _ExceptionHandler[_R_co]: ...

@type_check_only
class _SocketIOServerOptions(TypedDict, total=False):
    client_manager: Incomplete
    logger: Logger | bool
    json: Incomplete
    async_handlers: bool
    always_connect: bool

@type_check_only
class _EngineIOServerConfig(TypedDict, total=False):
    async_mode: Literal["threading", "eventlet", "gevent", "gevent_uwsgi"]
    ping_interval: float | tuple[float, float]  # seconds
    ping_timeout: float  # seconds
    max_http_buffer_size: int
    allow_upgrades: bool
    http_compression: bool
    compression_threshold: int
    cookie: str | dict[str, Any] | None
    cors_allowed_origins: str | list[str]
    cors_credentials: bool
    monitor_clients: bool
    engineio_logger: Logger | bool

@type_check_only
class _SocketIOKwargs(_SocketIOServerOptions, _EngineIOServerConfig): ...

class SocketIO:
    # This is an alias for `socketio.Server.reason` in `python-socketio`, which is not typed.
    reason: Incomplete
    # Many instance attributes are deliberately not included here,
    # as the maintainer of Flask-SocketIO considers them private, internal details:
    # https://github.com/python/typeshed/pull/10735#discussion_r1330768869
    def __init__(
        self,
        app: Flask | None = None,
        *,
        # SocketIO options
        manage_session: bool = True,
        message_queue: str | None = None,
        channel: str = "flask-socketio",
        path: str = "socket.io",
        resource: str = "socket.io",
        **kwargs: Unpack[_SocketIOKwargs],
    ) -> None: ...
    def init_app(
        self,
        app: Flask,
        *,
        # SocketIO options
        manage_session: bool = True,
        message_queue: str | None = None,
        channel: str = "flask-socketio",
        path: str = "socket.io",
        resource: str = "socket.io",
        **kwargs: Unpack[_SocketIOKwargs],
    ) -> None: ...
    def on(self, message: str, namespace: str | None = None) -> _HandlerDecorator: ...
    def on_error(self, namespace: str | None = None) -> _ExceptionHandlerDecorator: ...
    def on_error_default(self, exception_handler: _ExceptionHandler[_R_co]) -> _ExceptionHandler[_R_co]: ...
    def on_event(self, message: str, handler: _Handler[[Incomplete], object], namespace: str | None = None) -> None: ...
    @overload
    def event(self, event_handler: _Handler[_P, _R_co], /) -> _Handler[_P, _R_co]: ...
    @overload
    def event(self, namespace: str | None = None, *args, **kwargs) -> _HandlerDecorator: ...
    def on_namespace(self, namespace_handler: Namespace) -> None: ...
    def emit(
        self,
        event: str,
        *args,
        namespace: str = "/",  # / is the default (global) namespace
        to: str | None = None,
        include_self: bool = True,
        skip_sid: str | list[str] | None = None,
        callback: Callable[..., Incomplete] | None = None,
    ) -> None: ...
    def call(
        self,
        event: str,
        *args,
        namespace: str = "/",  # / is the default (global) namespace
        to: str | None = None,
        timeout: int = 60,  # seconds
        ignore_queue: bool = False,
    ): ...
    def send(
        self,
        data: Any,
        json: bool = False,
        namespace: str | None = None,
        to: str | None = None,
        callback: Callable[..., Incomplete] | None = None,
        include_self: bool = True,
        skip_sid: list[str] | str | None = None,
        **kwargs,
    ) -> None: ...
    def close_room(self, room: str, namespace: str | None = None) -> None: ...
    def run(
        self,
        app,
        host: str | None = None,
        port: int | None = None,
        *,
        debug: bool = True,
        use_reloader: bool = ...,
        reloader_options: dict[str, Incomplete] = {},
        log_output: bool = ...,
        allow_unsafe_werkzeug: bool = False,
        **kwargs,
    ) -> None: ...
    def stop(self) -> None: ...
    def start_background_task(self, target: Callable[_P, None], *args: _P.args, **kwargs: _P.kwargs) -> Thread: ...
    def sleep(self, seconds: int = 0): ...
    def test_client(
        self,
        app: Flask,
        namespace: str | None = None,
        query_string: str | None = None,
        headers: dict[str, Incomplete] | None = None,
        auth: dict[str, Incomplete] | None = None,
        flask_test_client: FlaskClient | None = None,
    ) -> SocketIOTestClient: ...

def emit(
    event,
    *args,
    namespace: str = "/",  # / is the default (global) namespace
    to: str | None = None,
    include_self: bool = True,
    skip_sid: str | list[str] | None = None,
    callback: Callable[..., Incomplete] | None = None,
    broadcast: bool = False,
) -> None: ...
def send(message: str, **kwargs) -> None: ...
def join_room(room, sid: str | None = None, namespace: str | None = None) -> None: ...
def leave_room(room, sid: str | None = None, namespace: str | None = None) -> None: ...
def close_room(room, namespace: str | None = None) -> None: ...
def rooms(sid: str | None = None, namespace: str | None = None) -> list[str]: ...
def disconnect(sid: str | None = None, namespace: str | None = None, silent: bool = False) -> None: ...
