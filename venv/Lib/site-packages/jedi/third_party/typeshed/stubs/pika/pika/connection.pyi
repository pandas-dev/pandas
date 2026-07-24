import abc
import ssl
from _typeshed import Incomplete
from collections.abc import Callable
from logging import Logger
from typing import Final
from typing_extensions import Self

from .callback import CallbackManager
from .channel import Channel
from .compat import AbstractBase
from .credentials import _Credentials
from .frame import Method
from .spec import Connection as SpecConnection

PRODUCT: str
LOGGER: Logger

class Parameters:
    __slots__ = (
        "_blocked_connection_timeout",
        "_channel_max",
        "_client_properties",
        "_connection_attempts",
        "_credentials",
        "_frame_max",
        "_heartbeat",
        "_host",
        "_locale",
        "_port",
        "_retry_delay",
        "_socket_timeout",
        "_stack_timeout",
        "_ssl_options",
        "_virtual_host",
        "_tcp_options",
    )
    DEFAULT_USERNAME: str
    DEFAULT_PASSWORD: str
    DEFAULT_BLOCKED_CONNECTION_TIMEOUT: None
    DEFAULT_CHANNEL_MAX: int
    DEFAULT_CLIENT_PROPERTIES: None
    DEFAULT_CREDENTIALS: Incomplete
    DEFAULT_CONNECTION_ATTEMPTS: int
    DEFAULT_FRAME_MAX: int
    DEFAULT_HEARTBEAT_TIMEOUT: None
    DEFAULT_HOST: str
    DEFAULT_LOCALE: str
    DEFAULT_PORT: int
    DEFAULT_RETRY_DELAY: float
    DEFAULT_SOCKET_TIMEOUT: float
    DEFAULT_STACK_TIMEOUT: float
    DEFAULT_SSL: bool
    DEFAULT_SSL_OPTIONS: None
    DEFAULT_SSL_PORT: int
    DEFAULT_VIRTUAL_HOST: str
    DEFAULT_TCP_OPTIONS: None
    def __init__(self) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def blocked_connection_timeout(self) -> float | None: ...
    @blocked_connection_timeout.setter
    def blocked_connection_timeout(self, value: float | None) -> None: ...
    @property
    def channel_max(self) -> int: ...
    @channel_max.setter
    def channel_max(self, value: int) -> None: ...
    @property
    def client_properties(self) -> dict[Incomplete, Incomplete] | None: ...
    @client_properties.setter
    def client_properties(self, value: dict[Incomplete, Incomplete] | None) -> None: ...
    @property
    def connection_attempts(self) -> int: ...
    @connection_attempts.setter
    def connection_attempts(self, value: int) -> None: ...
    @property
    def credentials(self) -> _Credentials: ...
    @credentials.setter
    def credentials(self, value: _Credentials) -> None: ...
    @property
    def frame_max(self) -> int: ...
    @frame_max.setter
    def frame_max(self, value: int) -> None: ...
    @property
    def heartbeat(self) -> int | Callable[[Connection, int], int] | None: ...
    @heartbeat.setter
    def heartbeat(self, value: int | Callable[[Connection, int], int] | None) -> None: ...
    @property
    def host(self) -> str: ...
    @host.setter
    def host(self, value: str) -> None: ...
    @property
    def locale(self) -> str: ...
    @locale.setter
    def locale(self, value: str) -> None: ...
    @property
    def port(self) -> int: ...
    @port.setter
    def port(self, value: int | str) -> None: ...
    @property
    def retry_delay(self) -> int | float: ...
    @retry_delay.setter
    def retry_delay(self, value: float) -> None: ...
    @property
    def socket_timeout(self) -> float | None: ...
    @socket_timeout.setter
    def socket_timeout(self, value: float | None) -> None: ...
    @property
    def stack_timeout(self) -> float | None: ...
    @stack_timeout.setter
    def stack_timeout(self, value: float | None) -> None: ...
    @property
    def ssl_options(self) -> SSLOptions | None: ...
    @ssl_options.setter
    def ssl_options(self, value: SSLOptions | None) -> None: ...
    @property
    def virtual_host(self) -> str: ...
    @virtual_host.setter
    def virtual_host(self, value: str) -> None: ...
    @property
    def tcp_options(self) -> dict[Incomplete, Incomplete] | None: ...
    @tcp_options.setter
    def tcp_options(self, value: dict[Incomplete, Incomplete] | None) -> None: ...

class ConnectionParameters(Parameters):
    __slots__ = ()
    def __init__(
        self,
        host: str = ...,
        port: int | str = ...,
        virtual_host: str = ...,
        credentials: _Credentials = ...,
        channel_max: int = ...,
        frame_max: int = ...,
        heartbeat: int | Callable[[Connection, int], int] | None = ...,
        ssl_options: SSLOptions | None = ...,
        connection_attempts: int = ...,
        retry_delay: float = ...,
        socket_timeout: float | None = ...,
        stack_timeout: float | None = ...,
        locale: str = ...,
        blocked_connection_timeout: float | None = ...,
        client_properties: dict[Incomplete, Incomplete] | None = ...,
        tcp_options: dict[Incomplete, Incomplete] | None = ...,
    ) -> None: ...

class URLParameters(Parameters):
    __slots__ = ("_all_url_query_values",)
    def __init__(self, url: str) -> None: ...

class SSLOptions:
    __slots__ = ("context", "server_hostname")
    context: ssl.SSLContext
    server_hostname: str | None
    def __init__(self, context: ssl.SSLContext, server_hostname: str | None = None) -> None: ...

class Connection(AbstractBase, metaclass=abc.ABCMeta):
    ON_CONNECTION_CLOSED: Final[str]
    ON_CONNECTION_ERROR: Final[str]
    ON_CONNECTION_OPEN_OK: Final[str]
    CONNECTION_CLOSED: Final[int]
    CONNECTION_INIT: Final[int]
    CONNECTION_PROTOCOL: Final[int]
    CONNECTION_START: Final[int]
    CONNECTION_TUNE: Final[int]
    CONNECTION_OPEN: Final[int]
    CONNECTION_CLOSING: Final[int]
    connection_state: int  # one of the constants above
    params: Parameters
    callbacks: CallbackManager
    server_capabilities: Incomplete
    server_properties: Incomplete
    known_hosts: Incomplete
    def __init__(
        self,
        parameters: Parameters | None = None,
        on_open_callback: Callable[[Self], object] | None = None,
        on_open_error_callback: Callable[[Self, BaseException], object] | None = None,
        on_close_callback: Callable[[Self, BaseException], object] | None = None,
        internal_connection_workflow: bool = True,
    ) -> None: ...
    def add_on_close_callback(self, callback: Callable[[Self, BaseException], object]) -> None: ...
    def add_on_connection_blocked_callback(self, callback: Callable[[Self, Method[SpecConnection.Blocked]], object]) -> None: ...
    def add_on_connection_unblocked_callback(
        self, callback: Callable[[Self, Method[SpecConnection.Unblocked]], object]
    ) -> None: ...
    def add_on_open_callback(self, callback: Callable[[Self], object]) -> None: ...
    def add_on_open_error_callback(
        self, callback: Callable[[Self, BaseException], object], remove_default: bool = True
    ) -> None: ...
    def channel(
        self, channel_number: int | None = None, on_open_callback: Callable[[Channel], object] | None = None
    ) -> Channel: ...
    def update_secret(self, new_secret, reason, callback=None) -> None: ...
    def close(self, reply_code: int = 200, reply_text: str = "Normal shutdown") -> None: ...
    @property
    def is_closed(self) -> bool: ...
    @property
    def is_closing(self) -> bool: ...
    @property
    def is_open(self) -> bool: ...
    @property
    def basic_nack(self) -> bool: ...
    @property
    def consumer_cancel_notify(self) -> bool: ...
    @property
    def exchange_exchange_bindings(self) -> bool: ...
    @property
    def publisher_confirms(self) -> bool: ...
