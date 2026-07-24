import datetime
from _typeshed import Incomplete, MaybeNone, Unused
from collections.abc import Callable, Iterator
from json import JSONDecoder
from typing import Any
from typing_extensions import Self, TypeAlias

from urllib3 import exceptions as urllib3_exceptions, fields, filepost, util
from urllib3.response import HTTPResponse

from . import auth, cookies, exceptions, hooks, status_codes, utils
from .adapters import HTTPAdapter
from .cookies import RequestsCookieJar
from .structures import CaseInsensitiveDict as CaseInsensitiveDict

_JSON: TypeAlias = Any  # any object that can be serialized to JSON

default_hooks = hooks.default_hooks
HTTPBasicAuth = auth.HTTPBasicAuth
cookiejar_from_dict = cookies.cookiejar_from_dict
get_cookie_header = cookies.get_cookie_header
RequestField = fields.RequestField
encode_multipart_formdata = filepost.encode_multipart_formdata
parse_url = util.parse_url
DecodeError = urllib3_exceptions.DecodeError
ReadTimeoutError = urllib3_exceptions.ReadTimeoutError
ProtocolError = urllib3_exceptions.ProtocolError
LocationParseError = urllib3_exceptions.LocationParseError
HTTPError = exceptions.HTTPError
MissingSchema = exceptions.MissingSchema
InvalidURL = exceptions.InvalidURL
ChunkedEncodingError = exceptions.ChunkedEncodingError
ContentDecodingError = exceptions.ContentDecodingError
ConnectionError = exceptions.ConnectionError
StreamConsumedError = exceptions.StreamConsumedError
guess_filename = utils.guess_filename
get_auth_from_url = utils.get_auth_from_url
requote_uri = utils.requote_uri
stream_decode_response_unicode = utils.stream_decode_response_unicode
to_key_val_list = utils.to_key_val_list
parse_header_links = utils.parse_header_links
iter_slices = utils.iter_slices
guess_json_utf = utils.guess_json_utf
super_len = utils.super_len
to_native_string = utils.to_native_string
codes = status_codes.codes

REDIRECT_STATI: Incomplete
DEFAULT_REDIRECT_LIMIT: Incomplete
CONTENT_CHUNK_SIZE: Incomplete
ITER_CHUNK_SIZE: Incomplete

class RequestEncodingMixin:
    @property
    def path_url(self) -> str: ...

class RequestHooksMixin:
    def register_hook(self, event, hook): ...
    def deregister_hook(self, event, hook): ...

class Request(RequestHooksMixin):
    hooks: Incomplete
    method: Incomplete
    url: Incomplete
    headers: Incomplete
    files: Incomplete
    data: Incomplete
    json: _JSON | None
    params: Incomplete
    auth: Incomplete
    cookies: Incomplete
    def __init__(
        self,
        method=None,
        url=None,
        headers=None,
        files=None,
        data=None,
        params=None,
        auth=None,
        cookies=None,
        hooks=None,
        json: _JSON | None = None,
    ) -> None: ...
    def prepare(self) -> PreparedRequest: ...

class PreparedRequest(RequestEncodingMixin, RequestHooksMixin):
    method: str | None
    url: str | None
    headers: CaseInsensitiveDict[str]
    body: bytes | str | None
    hooks: Incomplete
    def __init__(self) -> None: ...
    def prepare(
        self,
        method=None,
        url=None,
        headers=None,
        files=None,
        data=None,
        params=None,
        auth=None,
        cookies=None,
        hooks=None,
        json=None,
    ) -> None: ...
    def copy(self) -> PreparedRequest: ...
    def prepare_method(self, method) -> None: ...
    def prepare_url(self, url, params) -> None: ...
    def prepare_headers(self, headers) -> None: ...
    def prepare_body(self, data, files, json=None) -> None: ...
    def prepare_content_length(self, body: bytes | str | None) -> None: ...
    def prepare_auth(self, auth, url="") -> None: ...
    def prepare_cookies(self, cookies) -> None: ...
    def prepare_hooks(self, hooks) -> None: ...

class Response:
    __attrs__: Incomplete
    _content: bytes | None  # undocumented
    status_code: int
    headers: CaseInsensitiveDict[str]
    raw: HTTPResponse | MaybeNone
    url: str
    encoding: str | None
    history: list[Response]
    reason: str
    cookies: RequestsCookieJar
    elapsed: datetime.timedelta
    request: PreparedRequest
    connection: HTTPAdapter
    def __init__(self) -> None: ...
    def __bool__(self) -> bool: ...
    def __nonzero__(self) -> bool: ...
    def __iter__(self) -> Iterator[bytes]: ...
    def __enter__(self) -> Self: ...
    def __exit__(self, *args: Unused) -> None: ...
    @property
    def next(self) -> PreparedRequest | None: ...
    @property
    def ok(self) -> bool: ...
    @property
    def is_redirect(self) -> bool: ...
    @property
    def is_permanent_redirect(self) -> bool: ...
    @property
    def apparent_encoding(self) -> str: ...
    def iter_content(self, chunk_size: int | None = 1, decode_unicode: bool = False) -> Iterator[Incomplete]: ...
    def iter_lines(
        self, chunk_size: int | None = 512, decode_unicode: bool = False, delimiter: str | bytes | None = None
    ) -> Iterator[Incomplete]: ...
    @property
    def content(self) -> bytes | MaybeNone: ...
    @property
    def text(self) -> str: ...
    def json(
        self,
        *,
        cls: type[JSONDecoder] | None = ...,
        object_hook: Callable[[dict[Any, Any]], Any] | None = ...,
        parse_float: Callable[[str], Any] | None = ...,
        parse_int: Callable[[str], Any] | None = ...,
        parse_constant: Callable[[str], Any] | None = ...,
        object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = ...,
        **kwds: Any,
    ) -> Any: ...
    @property
    def links(self) -> dict[Incomplete, Incomplete]: ...
    def raise_for_status(self) -> None: ...
    def close(self) -> None: ...
