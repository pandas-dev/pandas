import email.message
import http.client
import re
from _ssl import _PasswordType
from _typeshed import Incomplete, MaybeNone, StrOrBytesPath
from collections.abc import Generator
from ssl import _SSLMethod
from typing import ClassVar, Final, Literal, TypeVar
from typing_extensions import Self

from .error import *

_T = TypeVar("_T", default=str)

__author__: Final[str]
__copyright__: Final[str]
__contributors__: Final[list[str]]
__license__: Final[str]
__version__: Final[str]

def has_timeout(timeout: float | None) -> bool: ...

debuglevel: Final[int]
RETRIES: Final[int]
DEFAULT_MAX_REDIRECTS: Final[int]
HOP_BY_HOP: Final[list[str]]
SAFE_METHODS: Final[tuple[str, ...]]
REDIRECT_CODES: Final[frozenset[int]]
CA_CERTS: Final[str]
DEFAULT_TLS_VERSION: Final[_SSLMethod]
URI: Final[re.Pattern[str]]

def parse_uri(uri: str) -> tuple[str | MaybeNone, str | MaybeNone, str | MaybeNone, str | MaybeNone, str | MaybeNone]: ...
def urlnorm(uri: str) -> tuple[str | MaybeNone, str | MaybeNone, str | MaybeNone, str | MaybeNone]: ...

re_url_scheme: Final[re.Pattern[str]]
re_unsafe: Final[re.Pattern[str]]

def safename(filename: str | bytes) -> str: ...

NORMALIZE_SPACE: Final[re.Pattern[str]]
USE_WWW_AUTH_STRICT_PARSING: Final[int]

class Authentication:
    path: Incomplete
    host: Incomplete
    credentials: Incomplete
    http: Incomplete
    def __init__(self, credentials, host, request_uri: str, headers, response, content, http) -> None: ...
    def depth(self, request_uri: str) -> int: ...
    def inscope(self, host: str, request_uri: str) -> bool: ...
    def request(self, method, request_uri, headers, content) -> None: ...
    def response(self, response, content) -> bool: ...
    def __eq__(self, auth: object) -> bool: ...
    def __ne__(self, auth: object) -> bool: ...
    def __lt__(self, auth: object) -> bool: ...
    def __gt__(self, auth: object) -> bool: ...
    def __le__(self, auth: object) -> bool: ...
    def __ge__(self, auth: object) -> bool: ...
    def __bool__(self) -> bool: ...

class BasicAuthentication(Authentication):
    def __init__(self, credentials, host, request_uri, headers, response, content, http) -> None: ...
    def request(self, method, request_uri, headers, content) -> None: ...

class DigestAuthentication(Authentication):
    challenge: Incomplete
    A1: Incomplete
    def __init__(self, credentials, host, request_uri, headers, response, content, http) -> None: ...
    def request(self, method, request_uri, headers, content, cnonce=None): ...
    def response(self, response, content) -> bool: ...

class HmacDigestAuthentication(Authentication):
    challenge: Incomplete
    hashmod: Incomplete
    pwhashmod: Incomplete
    key: Incomplete
    __author__: ClassVar[str]
    def __init__(self, credentials, host, request_uri, headers, response, content, http) -> None: ...
    def request(self, method, request_uri, headers, content) -> None: ...
    def response(self, response, content) -> bool: ...

class WsseAuthentication(Authentication):
    def __init__(self, credentials, host, request_uri, headers, response, content, http) -> None: ...
    def request(self, method, request_uri, headers, content) -> None: ...

class GoogleLoginAuthentication(Authentication):
    Auth: str
    def __init__(self, credentials, host, request_uri, headers, response, content, http) -> None: ...
    def request(self, method, request_uri, headers, content) -> None: ...

class FileCache:
    cache: Incomplete
    safe: Incomplete
    def __init__(self, cache, safe=...) -> None: ...
    def get(self, key): ...
    def set(self, key, value) -> None: ...
    def delete(self, key) -> None: ...

class Credentials:
    credentials: Incomplete
    def __init__(self) -> None: ...
    def add(self, name, password, domain: str = "") -> None: ...
    def clear(self) -> None: ...
    def iter(self, domain) -> Generator[tuple[str, str]]: ...

class KeyCerts(Credentials):
    def add(self, key, cert, domain, password) -> None: ...  # type: ignore[override]
    def iter(self, domain) -> Generator[tuple[str, str, str]]: ...  # type: ignore[override]

class AllHosts: ...

class ProxyInfo:
    bypass_hosts: Incomplete
    def __init__(
        self, proxy_type, proxy_host, proxy_port, proxy_rdns: bool = True, proxy_user=None, proxy_pass=None, proxy_headers=None
    ) -> None: ...
    def astuple(self): ...
    def isgood(self): ...
    def applies_to(self, hostname): ...
    def bypass_host(self, hostname): ...

def proxy_info_from_environment(method: Literal["http", "https"] = "http") -> ProxyInfo | None: ...
def proxy_info_from_url(url: str, method: Literal["http", "https"] = "http", noproxy: str | None = None) -> ProxyInfo: ...

class HTTPConnectionWithTimeout(http.client.HTTPConnection):
    proxy_info: Incomplete
    def __init__(self, host, port=None, timeout=None, proxy_info=None) -> None: ...
    sock: Incomplete
    def connect(self) -> None: ...

class HTTPSConnectionWithTimeout(http.client.HTTPSConnection):
    disable_ssl_certificate_validation: bool
    ca_certs: StrOrBytesPath | None
    proxy_info: Incomplete
    key_file: StrOrBytesPath | None
    cert_file: StrOrBytesPath | None
    key_password: _PasswordType | None
    def __init__(
        self,
        host: str,
        port: int | None = None,
        key_file: StrOrBytesPath | None = None,
        cert_file: StrOrBytesPath | None = None,
        timeout: float | None = None,
        proxy_info=None,
        ca_certs: StrOrBytesPath | None = None,
        disable_ssl_certificate_validation: bool = False,
        tls_maximum_version=None,
        tls_minimum_version=None,
        key_password: _PasswordType | None = None,
    ) -> None: ...
    sock: Incomplete
    def connect(self) -> None: ...

SCHEME_TO_CONNECTION: Final[dict[Literal["http", "https"], type[http.client.HTTPConnection]]]

class Http:
    proxy_info: Incomplete
    ca_certs: Incomplete
    disable_ssl_certificate_validation: bool
    tls_maximum_version: Incomplete
    tls_minimum_version: Incomplete
    connections: Incomplete
    cache: FileCache
    credentials: Credentials
    certificates: KeyCerts
    authorizations: list[Authentication]
    follow_redirects: bool
    redirect_codes: frozenset[int]
    optimistic_concurrency_methods: list[str]
    safe_methods: list[str]
    follow_all_redirects: bool
    ignore_etag: bool
    force_exception_to_status_code: bool
    timeout: float | None
    forward_authorization_headers: bool
    def __init__(
        self,
        cache: str | FileCache | None = None,
        timeout: float | None = None,
        proxy_info=...,
        ca_certs=None,
        disable_ssl_certificate_validation: bool = False,
        tls_maximum_version=None,
        tls_minimum_version=None,
    ) -> None: ...
    def close(self) -> None: ...
    def add_credentials(self, name, password, domain: str = "") -> None: ...
    def add_certificate(self, key, cert, domain, password=None) -> None: ...
    def clear_credentials(self) -> None: ...
    def request(self, uri, method: str = "GET", body=None, headers=None, redirections=5, connection_type=None): ...

class Response(dict[str, str | _T]):
    fromcache: bool
    version: int
    status: int
    reason: str
    previous: Response[_T] | None
    def __init__(self, info: http.client.HTTPResponse | email.message.Message | dict[str, _T]) -> None: ...
    @property
    def dict(self) -> Self: ...

__all__ = [
    "debuglevel",
    "FailedToDecompressContent",
    "Http",
    "HttpLib2Error",
    "ProxyInfo",
    "RedirectLimit",
    "RedirectMissingLocation",
    "Response",
    "RETRIES",
    "UnimplementedDigestAuthOptionError",
    "UnimplementedHmacDigestAuthOptionError",
]
