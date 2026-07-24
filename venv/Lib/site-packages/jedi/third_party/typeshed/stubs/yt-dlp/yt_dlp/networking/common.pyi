import abc
import enum
import io
from _typeshed import Unused
from collections.abc import Callable, Iterable, Mapping
from email.message import Message
from logging import Logger
from typing import IO, Any, Final
from typing_extensions import Self, TypeAlias

from ..cookies import YoutubeDLCookieJar
from ..utils._utils import _YDLLogger
from ..utils.networking import HTTPHeaderDict

DEFAULT_TIMEOUT: Final = 20
_RequestData: TypeAlias = bytes | Iterable[bytes] | IO[Any] | None
_Preference: TypeAlias = Callable[[RequestHandler, Request], int]

def register_preference(*handlers: type[RequestHandler]) -> Callable[..., object]: ...

class RequestDirector:
    handlers: dict[str, RequestHandler]
    preferences: set[_Preference]
    logger: Logger
    verbose: bool
    def __init__(self, logger: Logger, verbose: bool = False) -> None: ...
    def close(self) -> None: ...
    def add_handler(self, handler: RequestHandler) -> None: ...
    def send(self, request: Request) -> Response: ...

def register_rh(handler: RequestHandler) -> RequestHandler: ...

class Features(enum.Enum):
    ALL_PROXY = 1
    NO_PROXY = 2

class RequestHandler(abc.ABC, metaclass=abc.ABCMeta):
    headers: HTTPHeaderDict | dict[str, str]
    cookiejar: YoutubeDLCookieJar | None
    timeout: float | int
    proxies: Mapping[str, Any] | dict[str, Any]
    source_address: str | None
    verbose: bool
    prefer_system_certs: bool
    verify: bool
    legacy_ssl_support: bool
    def __init__(
        self,
        *,
        logger: _YDLLogger,
        headers: HTTPHeaderDict | Mapping[str, str] | None = None,
        cookiejar: YoutubeDLCookieJar | None = None,
        timeout: float | None = None,
        proxies: Mapping[str, Any] | None = None,
        source_address: str | None = None,
        verbose: bool = False,
        prefer_system_certs: bool = False,
        client_cert: dict[str, str | None] | None = None,
        verify: bool = True,
        legacy_ssl_support: bool = False,
        **_: Unused,
    ) -> None: ...
    def validate(self, request: Request) -> None: ...
    def send(self, request: Request) -> Response: ...
    def close(self) -> None: ...
    @property
    def RH_NAME(cls) -> str: ...
    @property
    def RH_KEY(cls) -> str: ...
    def __enter__(self) -> Self: ...
    def __exit__(self, *args: object) -> None: ...

class Request:
    proxies: Mapping[str, Any] | dict[str, Any]
    extensions: Mapping[str, Any] | dict[str, Any]
    def __init__(
        self,
        url: str,
        data: _RequestData | None = None,
        headers: HTTPHeaderDict | Mapping[str, str] | None = None,
        proxies: Mapping[str, Any] | None = None,
        query: Mapping[str, str] | None = None,
        method: str | None = None,
        extensions: Mapping[str, Any] | None = None,
    ) -> None: ...
    @property
    def url(self) -> str: ...
    @url.setter
    def url(self, url: str) -> None: ...
    @property
    def method(self) -> str: ...
    @method.setter
    def method(self, method: str) -> None: ...
    @property
    def data(self) -> _RequestData | io.IOBase: ...
    @data.setter
    def data(self, data: _RequestData) -> None: ...
    @property
    def headers(self) -> HTTPHeaderDict | dict[str, str]: ...
    @headers.setter
    def headers(self, new_headers: Mapping[str, str] | HTTPHeaderDict) -> None: ...
    def update(
        self,
        url: str | None = None,
        data: str | None = None,
        headers: HTTPHeaderDict | Mapping[str, str] | None = None,
        query: Mapping[str, str] | None = None,
        extensions: Mapping[str, Any] | None = None,
    ) -> None: ...
    def copy(self) -> Self: ...

def HEADRequest(
    url: str,
    data: _RequestData | None = None,
    headers: HTTPHeaderDict | Mapping[str, str] | None = None,
    proxies: Mapping[str, Any] | None = None,
    query: Mapping[str, str] | None = None,
    *,
    method: str = "HEAD",
    extensions: Mapping[str, Any] | None = None,
) -> Request: ...
def PATCHRequest(
    url: str,
    data: _RequestData | None = None,
    headers: HTTPHeaderDict | Mapping[str, str] | None = None,
    proxies: Mapping[str, Any] | None = None,
    query: Mapping[str, str] | None = None,
    *,
    method: str = "PATCH",
    extensions: Mapping[str, Any] | None = None,
) -> Request: ...
def PUTRequest(
    url: str,
    data: _RequestData | None = None,
    headers: HTTPHeaderDict | Mapping[str, str] | None = None,
    proxies: Mapping[str, Any] | None = None,
    query: Mapping[str, str] | None = None,
    *,
    method: str = "PUT",
    extensions: Mapping[str, Any] | None = None,
) -> Request: ...

class Response(io.IOBase):
    fp: io.IOBase
    headers: Message
    status: int
    url: str
    reason: str | None
    extensions: Mapping[str, Any] | dict[str, Any]
    def __init__(
        self,
        fp: io.IOBase,
        url: str,
        headers: Mapping[str, str],
        status: int = 200,
        reason: str | None = None,
        extensions: Mapping[str, Any] | dict[str, Any] | None = None,
    ) -> None: ...
    def readable(self) -> bool: ...
    def read(self, amt: int | None = None) -> bytes: ...
    def close(self) -> None: ...
    def get_header(self, name: str, default: str | None = None) -> str | None: ...
