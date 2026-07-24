from _typeshed import SupportsItems, SupportsRead, Unused
from collections.abc import Callable, Iterable, Mapping, MutableMapping
from typing import Any, TypedDict, type_check_only
from typing_extensions import Self, TypeAlias

from . import adapters, auth as _auth, compat, cookies, exceptions, hooks, models, status_codes, utils
from .models import _JSON, Response
from .structures import CaseInsensitiveDict as CaseInsensitiveDict

_BaseAdapter: TypeAlias = adapters.BaseAdapter
OrderedDict = compat.OrderedDict
cookiejar_from_dict = cookies.cookiejar_from_dict
extract_cookies_to_jar = cookies.extract_cookies_to_jar
RequestsCookieJar = cookies.RequestsCookieJar
merge_cookies = cookies.merge_cookies
Request = models.Request
PreparedRequest = models.PreparedRequest
DEFAULT_REDIRECT_LIMIT = models.DEFAULT_REDIRECT_LIMIT
default_hooks = hooks.default_hooks
dispatch_hook = hooks.dispatch_hook
to_key_val_list = utils.to_key_val_list
default_headers = utils.default_headers
to_native_string = utils.to_native_string
TooManyRedirects = exceptions.TooManyRedirects
InvalidSchema = exceptions.InvalidSchema
ChunkedEncodingError = exceptions.ChunkedEncodingError
ContentDecodingError = exceptions.ContentDecodingError
HTTPAdapter = adapters.HTTPAdapter
requote_uri = utils.requote_uri
get_environ_proxies = utils.get_environ_proxies
get_netrc_auth = utils.get_netrc_auth
should_bypass_proxies = utils.should_bypass_proxies
get_auth_from_url = utils.get_auth_from_url
codes = status_codes.codes
REDIRECT_STATI = models.REDIRECT_STATI

def preferred_clock() -> float: ...
def merge_setting(request_setting, session_setting, dict_class=...): ...
def merge_hooks(request_hooks, session_hooks, dict_class=...): ...

class SessionRedirectMixin:
    def resolve_redirects(
        self,
        resp,
        req,
        stream: bool = False,
        timeout=None,
        verify: bool = True,
        cert=None,
        proxies=None,
        yield_requests: bool = False,
        **adapter_kwargs,
    ): ...
    def rebuild_auth(self, prepared_request, response): ...
    def rebuild_proxies(self, prepared_request, proxies): ...
    def should_strip_auth(self, old_url, new_url): ...
    def rebuild_method(self, prepared_request: PreparedRequest, response: Response) -> None: ...
    def get_redirect_target(self, resp: Response) -> str | None: ...

_Data: TypeAlias = (
    # used in requests.models.PreparedRequest.prepare_body
    #
    # case: is_stream
    # see requests.adapters.HTTPAdapter.send
    # will be sent directly to http.HTTPConnection.send(...) (through urllib3)
    Iterable[bytes]
    # case: not is_stream
    # will be modified before being sent to urllib3.HTTPConnectionPool.urlopen(body=...)
    # see requests.models.RequestEncodingMixin._encode_params
    # see requests.models.RequestEncodingMixin._encode_files
    # note that keys&values are converted from Any to str by urllib.parse.urlencode
    | str
    | bytes
    | SupportsRead[str | bytes]
    | list[tuple[Any, Any]]
    | tuple[tuple[Any, Any], ...]
    | Mapping[Any, Any]
)
_Auth: TypeAlias = tuple[str, str] | _auth.AuthBase | Callable[[PreparedRequest], PreparedRequest]
_Cert: TypeAlias = str | tuple[str, str]
# Files is passed to requests.utils.to_key_val_list()
_FileName: TypeAlias = str | None
_FileContent: TypeAlias = SupportsRead[str | bytes] | str | bytes
_FileContentType: TypeAlias = str
_FileCustomHeaders: TypeAlias = Mapping[str, str]
_FileSpecTuple2: TypeAlias = tuple[_FileName, _FileContent]
_FileSpecTuple3: TypeAlias = tuple[_FileName, _FileContent, _FileContentType]
_FileSpecTuple4: TypeAlias = tuple[_FileName, _FileContent, _FileContentType, _FileCustomHeaders]
_FileSpec: TypeAlias = _FileContent | _FileSpecTuple2 | _FileSpecTuple3 | _FileSpecTuple4
_Files: TypeAlias = Mapping[str, _FileSpec] | Iterable[tuple[str, _FileSpec]]
_Hook: TypeAlias = Callable[[Response], Any]
_HooksInput: TypeAlias = Mapping[str, Iterable[_Hook] | _Hook]

_ParamsMappingKeyType: TypeAlias = str | bytes | int | float
_ParamsMappingValueType: TypeAlias = str | bytes | int | float | Iterable[str | bytes | int | float] | None
_Params: TypeAlias = (
    SupportsItems[_ParamsMappingKeyType, _ParamsMappingValueType]
    | tuple[_ParamsMappingKeyType, _ParamsMappingValueType]
    | Iterable[tuple[_ParamsMappingKeyType, _ParamsMappingValueType]]
    | str
    | bytes
)
_TextMapping: TypeAlias = MutableMapping[str, str]
_HeadersUpdateMapping: TypeAlias = Mapping[str, str | bytes | None]
_Timeout: TypeAlias = float | tuple[float | None, float | None]
_Verify: TypeAlias = bool | str

@type_check_only
class _Settings(TypedDict):
    verify: _Verify | None
    proxies: _TextMapping
    stream: bool
    cert: _Cert | None

class Session(SessionRedirectMixin):
    __attrs__: Any
    # See https://github.com/psf/requests/issues/5020#issuecomment-989082461:
    # requests sets this as a CaseInsensitiveDict, but users may set it to any MutableMapping
    headers: MutableMapping[str, str | bytes]
    auth: _Auth | None
    proxies: _TextMapping
    # Don't complain if:
    #   - value is assumed to be a list (which it is by default)
    #   - a _Hook is assigned directly, without wrapping it in a list (also works)
    hooks: dict[str, list[_Hook] | Any]
    params: _Params
    stream: bool
    verify: _Verify | None
    cert: _Cert | None
    max_redirects: int
    trust_env: bool
    cookies: RequestsCookieJar
    adapters: MutableMapping[str, adapters.BaseAdapter]
    def __init__(self) -> None: ...
    def __enter__(self) -> Self: ...
    def __exit__(self, *args: Unused) -> None: ...
    def prepare_request(self, request: Request) -> PreparedRequest: ...
    def request(
        self,
        method: str | bytes,
        url: str | bytes,
        params: _Params | None = None,
        data: _Data | None = None,
        headers: _HeadersUpdateMapping | None = None,
        cookies: None | RequestsCookieJar | _TextMapping = None,
        files: _Files | None = None,
        auth: _Auth | None = None,
        timeout: _Timeout | None = None,
        allow_redirects: bool = True,
        proxies: _TextMapping | None = None,
        hooks: _HooksInput | None = None,
        stream: bool | None = None,
        verify: _Verify | None = None,
        cert: _Cert | None = None,
        json: _JSON | None = None,
    ) -> Response: ...
    def get(
        self,
        url: str | bytes,
        *,
        params: _Params | None = ...,
        data: _Data | None = ...,
        headers: _HeadersUpdateMapping | None = ...,
        cookies: RequestsCookieJar | _TextMapping | None = ...,
        files: _Files | None = ...,
        auth: _Auth | None = ...,
        timeout: _Timeout | None = ...,
        allow_redirects: bool = ...,
        proxies: _TextMapping | None = ...,
        hooks: _HooksInput | None = ...,
        stream: bool | None = ...,
        verify: _Verify | None = ...,
        cert: _Cert | None = ...,
        json: _JSON | None = None,
    ) -> Response: ...
    def options(
        self,
        url: str | bytes,
        *,
        params: _Params | None = ...,
        data: _Data | None = ...,
        headers: _HeadersUpdateMapping | None = ...,
        cookies: RequestsCookieJar | _TextMapping | None = ...,
        files: _Files | None = ...,
        auth: _Auth | None = ...,
        timeout: _Timeout | None = ...,
        allow_redirects: bool = ...,
        proxies: _TextMapping | None = ...,
        hooks: _HooksInput | None = ...,
        stream: bool | None = ...,
        verify: _Verify | None = ...,
        cert: _Cert | None = ...,
        json: _JSON | None = None,
    ) -> Response: ...
    def head(
        self,
        url: str | bytes,
        *,
        params: _Params | None = ...,
        data: _Data | None = ...,
        headers: _HeadersUpdateMapping | None = ...,
        cookies: RequestsCookieJar | _TextMapping | None = ...,
        files: _Files | None = ...,
        auth: _Auth | None = ...,
        timeout: _Timeout | None = ...,
        allow_redirects: bool = ...,
        proxies: _TextMapping | None = ...,
        hooks: _HooksInput | None = ...,
        stream: bool | None = ...,
        verify: _Verify | None = ...,
        cert: _Cert | None = ...,
        json: _JSON | None = None,
    ) -> Response: ...
    def post(
        self,
        url: str | bytes,
        data: _Data | None = None,
        json: _JSON | None = None,
        *,
        params: _Params | None = ...,
        headers: _HeadersUpdateMapping | None = ...,
        cookies: RequestsCookieJar | _TextMapping | None = ...,
        files: _Files | None = ...,
        auth: _Auth | None = ...,
        timeout: _Timeout | None = ...,
        allow_redirects: bool = ...,
        proxies: _TextMapping | None = ...,
        hooks: _HooksInput | None = ...,
        stream: bool | None = ...,
        verify: _Verify | None = ...,
        cert: _Cert | None = ...,
    ) -> Response: ...
    def put(
        self,
        url: str | bytes,
        data: _Data | None = None,
        *,
        params: _Params | None = ...,
        headers: _HeadersUpdateMapping | None = ...,
        cookies: RequestsCookieJar | _TextMapping | None = ...,
        files: _Files | None = ...,
        auth: _Auth | None = ...,
        timeout: _Timeout | None = ...,
        allow_redirects: bool = ...,
        proxies: _TextMapping | None = ...,
        hooks: _HooksInput | None = ...,
        stream: bool | None = ...,
        verify: _Verify | None = ...,
        cert: _Cert | None = ...,
        json: _JSON | None = None,
    ) -> Response: ...
    def patch(
        self,
        url: str | bytes,
        data: _Data | None = None,
        *,
        params: _Params | None = ...,
        headers: _HeadersUpdateMapping | None = ...,
        cookies: RequestsCookieJar | _TextMapping | None = ...,
        files: _Files | None = ...,
        auth: _Auth | None = ...,
        timeout: _Timeout | None = ...,
        allow_redirects: bool = ...,
        proxies: _TextMapping | None = ...,
        hooks: _HooksInput | None = ...,
        stream: bool | None = ...,
        verify: _Verify | None = ...,
        cert: _Cert | None = ...,
        json: _JSON | None = None,
    ) -> Response: ...
    def delete(
        self,
        url: str | bytes,
        *,
        params: _Params | None = ...,
        data: _Data | None = ...,
        headers: _HeadersUpdateMapping | None = ...,
        cookies: RequestsCookieJar | _TextMapping | None = ...,
        files: _Files | None = ...,
        auth: _Auth | None = ...,
        timeout: _Timeout | None = ...,
        allow_redirects: bool = ...,
        proxies: _TextMapping | None = ...,
        hooks: _HooksInput | None = ...,
        stream: bool | None = ...,
        verify: _Verify | None = ...,
        cert: _Cert | None = ...,
        json: _JSON | None = None,
    ) -> Response: ...
    def send(
        self,
        request: PreparedRequest,
        *,
        stream: bool | None = ...,
        verify: _Verify | None = ...,
        proxies: _TextMapping | None = ...,
        cert: _Cert | None = ...,
        timeout: _Timeout | None = ...,
        allow_redirects: bool = ...,
        **kwargs: Any,
    ) -> Response: ...
    def merge_environment_settings(
        self,
        url: str | bytes | None,
        proxies: _TextMapping | None,
        stream: bool | None,
        verify: _Verify | None,
        cert: _Cert | None,
    ) -> _Settings: ...
    def get_adapter(self, url: str) -> _BaseAdapter: ...
    def close(self) -> None: ...
    def mount(self, prefix: str | bytes, adapter: _BaseAdapter) -> None: ...

def session() -> Session: ...
