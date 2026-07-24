import datetime
import io
from _typeshed import OptExcInfo, SupportsKeysAndGetItem, SupportsNoArgReadline, SupportsRead, WriteableBuffer
from _typeshed.wsgi import WSGIApplication, WSGIEnvironment
from collections.abc import Iterable, Mapping
from re import Pattern
from typing import IO, Any, ClassVar, Literal, Protocol, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import Self, TypeAlias

from webob._types import AsymmetricProperty, AsymmetricPropertyWithDelete, SymmetricProperty, SymmetricPropertyWithDelete
from webob.acceptparse import _AcceptCharsetProperty, _AcceptEncodingProperty, _AcceptLanguageProperty, _AcceptProperty
from webob.byterange import Range
from webob.cachecontrol import CacheControl
from webob.client import SendRequest
from webob.compat import cgi_FieldStorage
from webob.cookies import RequestCookies
from webob.descriptors import _authorization, _DateProperty
from webob.etag import IfRange, IfRangeDate, _ETagProperty
from webob.headers import EnvironHeaders
from webob.multidict import GetDict, MultiDict, NestedMultiDict, NoVars
from webob.response import Response

__all__ = ["BaseRequest", "Request", "LegacyRequest"]

_T = TypeVar("_T")
_HTTPMethod: TypeAlias = Literal["GET", "HEAD", "POST", "PUT", "DELETE", "CONNECT", "OPTIONS", "TRACE", "PATCH"]
_ListOrTuple: TypeAlias = list[_T] | tuple[_T, ...]
_RequestCacheControl: TypeAlias = CacheControl[Literal["request"]]

@type_check_only
class _SupportsReadAndNoArgReadline(SupportsRead[str | bytes], SupportsNoArgReadline[str | bytes], Protocol): ...

@type_check_only
class _RequestCacheControlDict(TypedDict, total=False):
    max_stale: int
    min_stale: int
    only_if_cached: bool
    no_cache: Literal[True] | str
    no_store: bool
    no_transform: bool
    max_age: int

_FieldStorageWithFile = cgi_FieldStorage

class _NoDefault: ...

NoDefault: _NoDefault

class BaseRequest:
    request_body_tempfile_limit: ClassVar[int]
    environ: WSGIEnvironment
    def __init__(self, environ: WSGIEnvironment, **kw: Any) -> None: ...
    @overload
    def encget(self, key: str, default: _T, encattr: str | None = None) -> str | _T: ...
    @overload
    def encget(self, key: str, *, encattr: str | None = None) -> str: ...
    def encset(self, key: str, val: str, encattr: str | None = None) -> None: ...
    @property
    def charset(self) -> str | None: ...
    def decode(self, charset: str | None = None, errors: str = "strict") -> Self: ...
    @property
    def body_file(self) -> SupportsRead[bytes]: ...
    @body_file.setter
    def body_file(self, value: SupportsRead[bytes]) -> None: ...
    @body_file.deleter
    def body_file(self) -> None: ...
    content_length: SymmetricPropertyWithDelete[int | None]
    body_file_raw: SymmetricProperty[SupportsRead[bytes]]
    is_body_seekable: bool
    @property
    def body_file_seekable(self) -> IO[bytes]: ...
    url_encoding: AsymmetricPropertyWithDelete[str, str | None]
    scheme: SymmetricProperty[str]
    method: AsymmetricPropertyWithDelete[_HTTPMethod, _HTTPMethod | None]
    http_version: SymmetricProperty[str]
    remote_user: SymmetricPropertyWithDelete[str | None]
    remote_host: SymmetricPropertyWithDelete[str | None]
    remote_addr: SymmetricPropertyWithDelete[str | None]
    query_string: AsymmetricPropertyWithDelete[str, str | None]
    server_name: SymmetricProperty[str]
    server_port: SymmetricProperty[int]
    script_name: AsymmetricPropertyWithDelete[str, str | None]
    path_info: SymmetricProperty[str]
    uscript_name = script_name  # bw compat
    upath_info = path_info  # bw compat
    content_type: AsymmetricPropertyWithDelete[str, str | None]
    headers: AsymmetricProperty[EnvironHeaders, SupportsKeysAndGetItem[str, str] | Iterable[tuple[str, str]]]
    @property
    def client_addr(self) -> str | None: ...
    @property
    def host_port(self) -> str: ...
    @property
    def host_url(self) -> str: ...
    @property
    def application_url(self) -> str: ...
    @property
    def path_url(self) -> str: ...
    @property
    def path(self) -> str: ...
    @property
    def path_qs(self) -> str: ...
    @property
    def url(self) -> str: ...
    def relative_url(self, other_url: str, to_application: bool = False) -> str: ...
    def path_info_pop(self, pattern: Pattern[str] | None = None) -> str | None: ...
    def path_info_peek(self) -> str | None: ...
    urlvars: SymmetricPropertyWithDelete[dict[str, str]]
    urlargs: SymmetricPropertyWithDelete[tuple[str, ...]]
    @property
    def is_xhr(self) -> bool: ...
    host: SymmetricPropertyWithDelete[str]
    @property
    def domain(self) -> str: ...
    @property
    def body(self) -> bytes: ...
    @body.setter
    def body(self, value: bytes | None) -> None: ...
    @body.deleter
    def body(self) -> None: ...
    json: SymmetricPropertyWithDelete[Any]
    json_body: SymmetricPropertyWithDelete[Any]
    text: SymmetricPropertyWithDelete[str]
    @property
    def POST(self) -> MultiDict[str, str | _FieldStorageWithFile] | NoVars: ...
    @property
    def GET(self) -> GetDict: ...
    @property
    def params(self) -> NestedMultiDict[str, str | _FieldStorageWithFile]: ...
    cookies: AsymmetricProperty[RequestCookies, SupportsKeysAndGetItem[str, str] | Iterable[tuple[str, str]]]
    def copy(self) -> Self: ...
    def copy_get(self) -> Self: ...
    @property
    def is_body_readable(self) -> bool: ...
    @is_body_readable.setter
    def is_body_readable(self, flag: bool) -> None: ...
    def make_body_seekable(self) -> None: ...
    def copy_body(self) -> None: ...
    def make_tempfile(self) -> io.BufferedRandom: ...
    def remove_conditional_headers(
        self, remove_encoding: bool = True, remove_range: bool = True, remove_match: bool = True, remove_modified: bool = True
    ) -> None: ...
    accept: _AcceptProperty
    accept_charset: _AcceptCharsetProperty
    accept_encoding: _AcceptEncodingProperty
    accept_language: _AcceptLanguageProperty
    authorization: AsymmetricPropertyWithDelete[_authorization | None, tuple[str, str | dict[str, str]] | list[Any] | str | None]
    cache_control: AsymmetricPropertyWithDelete[
        _RequestCacheControl, _RequestCacheControl | _RequestCacheControlDict | str | None
    ]
    if_match: _ETagProperty
    if_none_match: _ETagProperty
    date: _DateProperty
    if_modified_since: _DateProperty
    if_unmodified_since: _DateProperty
    if_range: AsymmetricPropertyWithDelete[
        IfRange | IfRangeDate, IfRange | IfRangeDate | datetime.datetime | datetime.date | str | None
    ]
    max_forwards: SymmetricPropertyWithDelete[int | None]
    pragma: SymmetricPropertyWithDelete[str | None]
    range: AsymmetricPropertyWithDelete[Range | None, tuple[int, int | None] | list[int | None] | list[int] | str | None]
    referer: SymmetricPropertyWithDelete[str | None]
    referrer = referer
    user_agent: SymmetricPropertyWithDelete[str | None]
    def as_bytes(self, skip_body: bool = False) -> bytes: ...
    def as_text(self) -> str: ...
    @classmethod
    def from_bytes(cls, b: bytes) -> Self: ...
    @classmethod
    def from_text(cls, s: str) -> Self: ...
    @classmethod
    def from_file(cls, fp: _SupportsReadAndNoArgReadline) -> Self: ...
    @overload
    def call_application(
        self, application: WSGIApplication, catch_exc_info: Literal[False] = False
    ) -> tuple[str, list[tuple[str, str]], Iterable[bytes]]: ...
    @overload
    def call_application(
        self, application: WSGIApplication, catch_exc_info: Literal[True]
    ) -> tuple[str, list[tuple[str, str]], Iterable[bytes], OptExcInfo | None]: ...
    @overload
    def call_application(
        self, application: WSGIApplication, catch_exc_info: bool
    ) -> (
        tuple[str, list[tuple[str, str]], Iterable[bytes], OptExcInfo | None] | tuple[str, list[tuple[str, str]], Iterable[bytes]]
    ): ...
    ResponseClass: type[Response]
    def send(self, application: WSGIApplication | None = None, catch_exc_info: bool = False) -> Response: ...
    get_response = send
    def make_default_send_app(self) -> SendRequest: ...
    @classmethod
    def blank(
        cls,
        path: str,
        environ: dict[str, None] | None = None,
        base_url: str | None = None,
        headers: Mapping[str, str] | None = None,
        POST: str | bytes | Mapping[Any, Any] | Mapping[Any, _ListOrTuple[Any]] | None = None,
        **kw: Any,
    ) -> Self: ...

class LegacyRequest(BaseRequest):
    @property  # type: ignore[override]
    def uscript_name(self) -> str: ...
    @uscript_name.setter
    def uscript_name(self, value: str) -> None: ...
    @property  # type: ignore[override]
    def upath_info(self) -> str: ...
    @upath_info.setter
    def upath_info(self, value: str) -> None: ...
    def encget(self, key: str, default: Any = ..., encattr: str | None = None) -> Any: ...

class AdhocAttrMixin:
    def __setattr__(self, attr: str, value: Any) -> None: ...
    def __getattr__(self, attr: str) -> Any: ...
    def __delattr__(self, attr: str) -> None: ...

class Request(AdhocAttrMixin, BaseRequest):
    # this is so Request doesn't count as callable, it's not very pretty
    # but we run into trouble with overlapping overloads in wsgify if we
    # don't exclude __call__ from arbitrary attribute access
    __call__: None

class DisconnectionError(IOError): ...

def environ_from_url(path: str) -> WSGIEnvironment: ...
def environ_add_POST(
    env: WSGIEnvironment,
    data: str | bytes | Mapping[Any, Any] | Mapping[Any, _ListOrTuple[Any]] | None,
    content_type: str | None = None,
) -> None: ...

class LimitedLengthFile(io.RawIOBase):
    file: SupportsRead[bytes]
    maxlen: int
    remaining: int
    def __init__(self, file: SupportsRead[bytes], maxlen: int) -> None: ...
    def fileno(self) -> int: ...
    @staticmethod
    def readable() -> Literal[True]: ...
    def readinto(self, buff: WriteableBuffer) -> int: ...

class Transcoder:
    charset: str
    errors: str
    def __init__(self, charset: str, errors: str = "strict") -> None: ...
    def transcode_query(self, q: str) -> str: ...
    def transcode_fs(self, fs: cgi_FieldStorage, content_type: str) -> io.BytesIO: ...
