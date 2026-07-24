from _typeshed import SupportsItems, SupportsRead
from _typeshed.wsgi import StartResponse, WSGIApplication, WSGIEnvironment
from collections.abc import Iterable, Iterator, Sequence
from datetime import timedelta
from typing import IO, Any, Literal, Protocol, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import Self, TypeAlias

from webob._types import AsymmetricProperty, AsymmetricPropertyWithDelete, SymmetricProperty, SymmetricPropertyWithDelete
from webob.byterange import ContentRange
from webob.cachecontrol import CacheControl
from webob.cookies import _SameSitePolicy
from webob.descriptors import _authorization, _ContentRangeParams, _DateProperty, _ListProperty
from webob.headers import ResponseHeaders
from webob.request import Request

__all__ = ["Response"]

_ResponseT = TypeVar("_ResponseT", bound=Response)
_ResponseCacheControl: TypeAlias = CacheControl[Literal["response"]]

@type_check_only
class _ResponseCacheExpires(Protocol):
    def __call__(
        self,
        seconds: int | timedelta = 0,
        *,
        public: bool = ...,
        private: Literal[True] | str = ...,
        no_cache: Literal[True] | str = ...,
        no_store: bool = ...,
        no_transform: bool = ...,
        must_revalidate: bool = ...,
        proxy_revalidate: bool = ...,
        max_age: int = ...,
        s_maxage: int = ...,
        s_max_age: int = ...,
        stale_while_revalidate: int = ...,
        stale_if_error: int = ...,
    ) -> None: ...

@type_check_only
class _ResponseCacheControlDict(TypedDict, total=False):
    public: bool
    private: Literal[True] | str
    no_cache: Literal[True] | str
    no_store: bool
    no_transform: bool
    must_revalidate: bool
    proxy_revalidate: bool
    max_age: int
    s_maxage: int
    s_max_age: int
    stale_while_revalidate: int
    stale_if_error: int

class Response:
    default_content_type: str
    default_charset: str
    unicode_errors: str
    default_conditional_response: bool
    default_body_encoding: str
    request: Request | None
    environ: WSGIEnvironment | None
    status: AsymmetricProperty[str, int | str | bytes]
    conditional_response: bool
    def __init__(
        self,
        body: bytes | str | None = None,
        status: int | str | bytes | None = None,
        headerlist: list[tuple[str, str]] | None = None,
        app_iter: Iterable[bytes] | None = None,
        content_type: str | None = None,
        conditional_response: bool | None = None,
        charset: str = ...,
        **kw: Any,
    ) -> None: ...
    @classmethod
    def from_file(cls, fp: IO[str] | IO[bytes]) -> Response: ...
    def copy(self) -> Response: ...
    status_code: SymmetricProperty[int]
    status_int: SymmetricProperty[int]
    headerlist: AsymmetricPropertyWithDelete[list[tuple[str, str]], Iterable[tuple[str, str]] | SupportsItems[str, str]]
    headers: AsymmetricProperty[ResponseHeaders, SupportsItems[str, str] | Iterable[tuple[str, str]]]
    body: SymmetricPropertyWithDelete[bytes]
    json: SymmetricPropertyWithDelete[Any]
    json_body: SymmetricPropertyWithDelete[Any]
    @property
    def has_body(self) -> bool: ...
    text: SymmetricPropertyWithDelete[str]
    unicode_body: SymmetricPropertyWithDelete[str]  # deprecated
    ubody: SymmetricPropertyWithDelete[str]  # deprecated
    body_file: AsymmetricPropertyWithDelete[ResponseBodyFile, SupportsRead[bytes]]
    content_length: AsymmetricPropertyWithDelete[int | None, int | str | bytes | None]
    def write(self, text: str | bytes) -> int: ...
    app_iter: SymmetricPropertyWithDelete[Iterable[bytes]]
    allow: _ListProperty
    vary: _ListProperty
    content_encoding: SymmetricPropertyWithDelete[str | None]
    content_language: SymmetricPropertyWithDelete[str | None]
    content_location: SymmetricPropertyWithDelete[str | None]
    content_md5: SymmetricPropertyWithDelete[str | None]
    content_disposition: SymmetricPropertyWithDelete[str | None]
    accept_ranges: SymmetricPropertyWithDelete[str | None]
    content_range: AsymmetricPropertyWithDelete[ContentRange | None, _ContentRangeParams]
    date: _DateProperty
    expires: _DateProperty
    last_modified: _DateProperty
    etag: AsymmetricPropertyWithDelete[str | None, tuple[str, bool] | str | None]
    @property
    def etag_strong(self) -> str | None: ...
    location: SymmetricPropertyWithDelete[str | None]
    pragma: SymmetricPropertyWithDelete[str | None]
    age: SymmetricPropertyWithDelete[int | None]
    retry_after: _DateProperty
    server: SymmetricPropertyWithDelete[str | None]
    www_authenticate: AsymmetricPropertyWithDelete[
        _authorization | None, tuple[str, str | dict[str, str]] | list[Any] | str | None
    ]
    charset: SymmetricPropertyWithDelete[str | None]
    content_type: SymmetricPropertyWithDelete[str | None]
    content_type_params: AsymmetricPropertyWithDelete[dict[str, str], SupportsItems[str, str] | None]
    def set_cookie(
        self,
        name: str | bytes,
        value: str | bytes | None = "",
        max_age: int | timedelta | None = None,
        path: str = "/",
        domain: str | None = None,
        secure: bool = False,
        httponly: bool = False,
        comment: str | None = None,
        overwrite: bool = False,
        samesite: _SameSitePolicy | None = None,
    ) -> None: ...
    def delete_cookie(self, name: str | bytes, path: str = "/", domain: str | None = None) -> None: ...
    def unset_cookie(self, name: str | bytes, strict: bool = True) -> None: ...
    @overload
    def merge_cookies(self, resp: _ResponseT) -> _ResponseT: ...
    @overload
    def merge_cookies(self, resp: WSGIApplication) -> WSGIApplication: ...
    cache_control: AsymmetricProperty[_ResponseCacheControl, _ResponseCacheControl | _ResponseCacheControlDict | str | None]
    cache_expires: AsymmetricProperty[_ResponseCacheExpires, timedelta | int | bool | None]
    def encode_content(self, encoding: Literal["gzip", "identity"] = "gzip", lazy: bool = False) -> None: ...
    def decode_content(self) -> None: ...
    def md5_etag(self, body: bytes | None = None, set_content_md5: bool = False) -> None: ...
    def __call__(self, environ: WSGIEnvironment, start_response: StartResponse) -> Iterable[bytes]: ...
    def conditional_response_app(self, environ: WSGIEnvironment, start_response: StartResponse) -> Iterable[bytes]: ...
    def app_iter_range(self, start: int, stop: int | None) -> AppIterRange: ...
    def __str__(self, skip_body: bool = False) -> str: ...

class ResponseBodyFile:
    mode: Literal["wb"]
    closed: Literal[False]
    response: Response
    def __init__(self, response: Response) -> None: ...
    @property
    def encoding(self) -> str | None: ...
    # NOTE: Technically this is an instance attribute and not a method
    def write(self, text: str | bytes) -> int: ...
    def writelines(self, seq: Sequence[str | bytes]) -> None: ...
    def flush(self) -> None: ...
    def tell(self) -> int: ...

class AppIterRange:
    app_iter: Iterator[bytes]
    start: int
    stop: int | None
    def __init__(self, app_iter: Iterable[bytes], start: int, stop: int | None) -> None: ...
    def __iter__(self) -> Self: ...
    def next(self) -> bytes: ...
    __next__ = next
    def close(self) -> None: ...

class EmptyResponse:
    def __init__(self, app_iter: Iterable[bytes] | None = None) -> None: ...
    def __iter__(self) -> Self: ...
    def __len__(self) -> Literal[0]: ...
    def next(self) -> bytes: ...
    __next__ = next
